#!/usr/bin/env python3

import argparse
import bisect
import fnmatch
import os
import shutil
import time
from pathlib import Path
from stat import S_ISDIR, S_ISLNK, S_ISREG

parser = argparse.ArgumentParser(
    description="Performs garbage collection for Snakemake file cache and pipeline scratch directories",
)

parser.add_argument("--snakemake-output-cache-dir", type=str, help="Will try to use $SNAKEMAKE_OUTPUT_CACHE by default")
parser.add_argument("--symlink-roots-dir", type=str, help="Cache entries that a symbolic link under this directory points to are alive and never removed")
parser.add_argument("--target-size", type=str, help="Target size of the cache after pruning (can use K,M,G suffixes)", required=True)
parser.add_argument("--min-free", type=str, help="After reaching the target size, keep removing entries (and inactive scratch directories, see --scratch-grace-period) until the filesystem holding the cache has at least this much free space (can use K,M,G suffixes)")
parser.add_argument("--grace-period", type=str, default="6h", help="Never remove entries modified more recently than this (can use s,m,h,d suffixes)")
parser.add_argument("--scratch-pattern", action="append", default=[], help="Glob for pipeline scratch directories directly under --symlink-roots-dir (repeatable), e.g. 'detector_benchmarks_*'")
parser.add_argument("--scratch-max-age", type=str, help="Always remove scratch directories inactive for longer than this (can use s,m,h,d suffixes)")
parser.add_argument("--scratch-grace-period", type=str, help="With --min-free, remove scratch directories inactive for longer than this, least recently active first, when removing cache entries alone does not free enough space (can use s,m,h,d suffixes)")

args = parser.parse_args()

if args.snakemake_output_cache_dir is None:
    if "SNAKEMAKE_OUTPUT_CACHE" in os.environ:
        cache_dir = Path(os.environ["SNAKEMAKE_OUTPUT_CACHE"])
    else:
        raise Exception("Must specify --snakemake-output-cache-dir or set $SNAKEMAKE_OUTPUT_CACHE")
else:
    cache_dir = Path(args.snakemake_output_cache_dir)
assert cache_dir.is_dir()
# symbolic links are compared after resolving them
cache_dir_resolved = Path(os.path.realpath(cache_dir))

def parse_duration(duration):
    if duration is None:
        return None
    if duration[-1] in ["s", "m", "h", "d"]:
        return float(duration[:-1]) * {
            "s": 1,
            "m": 60,
            "h": 3600,
            "d": 86400,
        }[duration[-1]]
    return float(duration)

def parse_size(size):
    if size[-1] in ["K", "M", "G"]:
        return int(size[:-1]) * {
            "K": 1024,
            "M": 1024**2,
            "G": 1024**3,
        }[size[-1]]
    else:
        return int(size)

grace_period = parse_duration(args.grace_period)
scratch_max_age = parse_duration(args.scratch_max_age)
scratch_grace_period = parse_duration(args.scratch_grace_period)
target_size = parse_size(args.target_size)
min_free = parse_size(args.min_free) if args.min_free is not None else None

if (scratch_max_age is not None or scratch_grace_period is not None) and (args.symlink_roots_dir is None or not args.scratch_pattern):
    parser.error("--scratch-max-age/--scratch-grace-period require --symlink-roots-dir and at least one --scratch-pattern")

def free_space():
    # space available to unprivileged users (statvfs f_bavail) on the filesystem holding the cache
    return shutil.disk_usage(cache_dir).free

GiB = 1024**3

def scan_tree(root):
    """Walk root without following symbolic links. Returns the newest mtime of
    root and anything below it, the allocated size, and the resolved targets of
    the symbolic links below it."""
    newest = 0.0
    allocated = 0
    targets = []
    stack = [root]
    while stack:
        path = stack.pop()
        try:
            st = os.lstat(path)
        except OSError:
            # removed meanwhile, or not accessible
            continue
        newest = max(newest, st.st_mtime)
        allocated += st.st_blocks * 512
        if S_ISLNK(st.st_mode):
            targets.append(Path(os.path.realpath(path)))
        elif S_ISDIR(st.st_mode):
            try:
                with os.scandir(path) as it:
                    stack.extend(entry.path for entry in it)
            except OSError:
                pass
    return newest, allocated, targets

def remove_tree(path):
    # rmtree never follows symbolic links; keep going on errors so that a
    # single unremovable file does not keep the rest of the directory
    def report(func, p, exc):
        print(f"  could not remove {p}: {exc[1] if isinstance(exc, tuple) else exc}")
    try:
        shutil.rmtree(path, onexc=report)
    except TypeError:
        # Python < 3.12
        shutil.rmtree(path, onerror=report)

# List the cache before scanning for symlinks: entries stored (and symlinked)
# while the scan is running are then never considered for removal
cache_paths = list(cache_dir.iterdir())

now = time.time()

# Scan the symlink roots. Symbolic links outside the pipeline scratch
# directories pin cache entries for good; those inside a scratch directory
# pin them only as long as that directory exists.
pinned_forever = set()
pins = {}  # resolved cache entry -> set of scratch directories with a symbolic link to it
scratch = []  # (path, last activity, allocated size, targets)
if args.symlink_roots_dir is not None:
    roots_dir = Path(args.symlink_roots_dir)
    for child in sorted(roots_dir.iterdir()):
        try:
            st = child.lstat()
        except OSError as e:
            # removed meanwhile, or not accessible: neither scanned nor removed
            if not isinstance(e, FileNotFoundError):
                print(f"Skipping {child}: {e}")
            continue
        is_scratch = (
            S_ISDIR(st.st_mode)
            and Path(os.path.realpath(child)) != cache_dir_resolved
            and any(fnmatch.fnmatch(child.name, p) for p in args.scratch_pattern)
        )
        last_activity, allocated, targets = scan_tree(child)
        if not is_scratch:
            pinned_forever.update(targets)
            continue
        inactive = now - last_activity
        if scratch_max_age is not None and inactive > scratch_max_age:
            print(f"Removing scratch directory {child} of {allocated / GiB:.1f} GiB, inactive for {inactive / 3600:.1f} h (max age {args.scratch_max_age})")
            remove_tree(child)
            continue
        scratch.append((child, last_activity, allocated, targets))
        for target in set(targets):
            pins.setdefault(target, set()).add(child)

def is_alive(key):
    return key in pinned_forever or bool(pins.get(key))

paths = []
total_size = 0
# space allocated on disk (st_blocks), unlike the apparent size (st_size) used for the target size
total_allocated = 0
grace_size = grace_count = 0
for cache_path in cache_paths:
    try:
        stat = cache_path.stat()
    except FileNotFoundError:
        # e.g. a snakemake temporary directory that has since been removed
        continue
    size = stat.st_size
    if S_ISREG(stat.st_mode) and size == 0:
        # empty entries are never valid (e.g. written on a full disk): remove regardless of age or symlinks
        cache_path.unlink(missing_ok=True)
        print(f"Removing empty {cache_path}")
        continue
    allocated = stat.st_blocks * 512
    total_size += size
    total_allocated += allocated
    if now - stat.st_mtime < grace_period:
        # recently stored, its symlinks may not have been seen by the scan
        grace_size += size
        grace_count += 1
        continue
    age = now - stat.st_atime
    score = size * age
    paths.append((score, cache_path, size, allocated))

# Removable entries sorted by increasing score (the next one to remove is
# last); alive entries are set aside until their last scratch directory is
# removed.
paths.sort(key=lambda t: t[0])
removable = []
waiting = {}  # resolved cache entry -> entry tuple
for entry in paths:
    key = cache_dir_resolved / entry[1].name
    if is_alive(key):
        waiting[key] = entry
    else:
        removable.append(entry)

print(f"Total cache size: {total_size / 1024. / 1024.:.1f} MiB")
if min_free is not None:
    print(f"Free space: {free_space() / GiB:.1f} GiB (minimum {min_free / GiB:.1f} GiB)")

def remove_next_entry():
    global total_size, total_allocated
    _, cache_path, size, allocated = removable.pop()
    cache_path.unlink(missing_ok=True)
    print(f"Removing {cache_path} of {size / 1024. / 1024.:.1f} MiB")
    total_size -= size
    total_allocated -= allocated

while total_size >= target_size and removable:
    remove_next_entry()

if min_free is not None:
    # The cache shares its filesystem with the pipeline scratch directories, so
    # a cache below the target size does not guarantee free space. Keep
    # removing entries in the same order. Free space is measured again after
    # every removal, because an entry that is still open (e.g. being read by a
    # running job) frees no space until it is closed. The grace period and the
    # alive protection still apply: removing those entries would break
    # running pipelines.
    #
    # When no removable entry is left, remove the least recently active
    # scratch directory that has been inactive for longer than the scratch
    # grace period (typically a failed pipeline that was never retried). The
    # cache entries that only it pointed to become removable.
    candidates = []
    if scratch_grace_period is not None:
        candidates = sorted(
            (s for s in scratch if now - s[1] > scratch_grace_period),
            key=lambda s: s[1],
            reverse=True,
        )
    while free_space() < min_free:
        if removable:
            remove_next_entry()
        elif candidates:
            path, last_activity, allocated, targets = candidates.pop()
            print(f"Removing scratch directory {path} of {allocated / GiB:.1f} GiB, inactive for {(now - last_activity) / 3600:.1f} h (scratch grace period {args.scratch_grace_period})")
            remove_tree(path)
            scratch.remove((path, last_activity, allocated, targets))
            for target in set(targets):
                pins[target].discard(path)
                if not is_alive(target) and target in waiting:
                    # tuples sort by score first, like removable
                    bisect.insort(removable, waiting.pop(target))
        else:
            break

    free = free_space()
    print(f"Total cache size after pruning: {total_size / 1024. / 1024.:.1f} MiB")
    print(f"Free space after pruning: {free / GiB:.1f} GiB (minimum {min_free / GiB:.1f} GiB)")
    if free < min_free:
        usage = shutil.disk_usage(cache_dir)
        alive_size = sum(e[2] for e in waiting.values())
        scratch_size = sum(s[2] for s in scratch)
        if scratch_grace_period is not None:
            scratch_note = f"{len(scratch)} scratch directories of {scratch_size / GiB:.1f} GiB are kept, all active less than {args.scratch_grace_period} ago (scratch grace period). "
        else:
            scratch_note = f"{len(scratch)} scratch directories of {scratch_size / GiB:.1f} GiB are kept (no --scratch-grace-period). "
        print(
            f"WARNING: free space on the filesystem holding {cache_dir} is {free / GiB:.1f} GiB, "
            f"below the minimum of {min_free / GiB:.1f} GiB (short by {(min_free - free) / GiB:.1f} GiB), "
            f"and no removable cache entries or scratch directories are left. "
            f"The remaining cache of {total_size / GiB:.1f} GiB is protected: "
            f"{grace_size / GiB:.1f} GiB in {grace_count} entries modified less than {args.grace_period} ago (grace period), "
            f"{alive_size / GiB:.1f} GiB in {len(waiting)} entries that symbolic links under {args.symlink_roots_dir or '(no --symlink-roots-dir)'} point to (alive). "
            f"{scratch_note}"
            f"Outside the cache, about {max(usage.used - total_allocated, 0) / GiB:.1f} GiB of the {usage.total / GiB:.1f} GiB filesystem is used "
            f"(scratch directories, other files, and removed cache entries or scratch files that are still open)."
        )
