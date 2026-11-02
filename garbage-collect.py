#!/usr/bin/env python3

import argparse
import os
import time
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Performs garbage collection for Snakemake file cache",
)

parser.add_argument("--snakemake-output-cache-dir", type=str, help="Will try to use $SNAKEMAKE_OUTPUT_CACHE by default")
parser.add_argument("--symlink-roots-dir", type=str)
parser.add_argument("--target-size", type=str, help="Target size of the cache after pruning (can use K,M,G suffixes)", required=True)
parser.add_argument("--grace-period", type=str, default="6h", help="Never remove entries modified more recently than this (can use s,m,h,d suffixes)")

args = parser.parse_args()

if args.snakemake_output_cache_dir is None:
    if "SNAKEMAKE_OUTPUT_CACHE" in os.environ:
        cache_dir = Path(os.environ["SNAKEMAKE_OUTPUT_CACHE"])
    else:
        raise Exception("Must specify --snakemake-output-cache-dir or set $SNAKEMAKE_OUTPUT_CACHE")
else:
    cache_dir = Path(args.snakemake_output_cache_dir)
assert cache_dir.is_dir()

if args.grace_period[-1] in ["s", "m", "h", "d"]:
    grace_period = float(args.grace_period[:-1]) * {
        "s": 1,
        "m": 60,
        "h": 3600,
        "d": 86400,
    }[args.grace_period[-1]]
else:
    grace_period = float(args.grace_period)

if args.target_size[-1] in ["K", "M", "G"]:
    target_size = int(args.target_size[:-1]) * {
        "K": 1024,
        "M": 1024**2,
        "G": 1024**3,
    }[args.target_size[-1]]
else:
    target_size = int(args.target_size)

# List the cache before scanning for symlinks: entries stored (and symlinked)
# while the scan is running are then never considered for removal
cache_paths = list(cache_dir.iterdir())

alive = set()
if args.symlink_roots_dir is not None:
    for path in Path(args.symlink_roots_dir).rglob("*"):
        if path.is_symlink():
            alive.add(path.resolve())

paths = []
total_size = 0
for cache_path in cache_paths:
    try:
        stat = cache_path.stat()
    except FileNotFoundError:
        # e.g. a snakemake temporary directory that has since been removed
        continue
    size = stat.st_size
    total_size += size
    if time.time() - stat.st_mtime < grace_period:
        # recently stored, its symlinks may not have been seen by the scan
        continue
    age = time.time() - stat.st_atime
    score = size * age
    paths.append((cache_path, score, size))

print(f"Total cache size: {total_size / 1024. / 1024.:.1f} MiB")

paths.sort(key=lambda t: t[1])
while total_size >= target_size and paths:
    cache_path, _, size = paths.pop()
    if cache_path in alive:
        print(f"{cache_path} is alive")
    else:
        cache_path.unlink(missing_ok=True)
        print(f"Removing {cache_path} of {size / 1024. / 1024.:.1f} MiB")
        total_size -= size
