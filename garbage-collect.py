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

args = parser.parse_args()

if args.snakemake_output_cache_dir is None:
    if "SNAKEMAKE_OUTPUT_CACHE" in os.environ:
        cache_dir = Path(os.environ["SNAKEMAKE_OUTPUT_CACHE"])
    else:
        raise Exception("Must specify --snakemake-output-cache-dir or set $SNAKEMAKE_OUTPUT_CACHE")
else:
    cache_dir = Path(args.snakemake_output_cache_dir)
assert cache_dir.is_dir()

if args.target_size[-1] in ["K", "M", "G"]:
    target_size = int(args.target_size[:-1]) * {
        "K": 1024,
        "M": 1024**2,
        "G": 1024**3,
    }[args.target_size[-1]]
else:
    target_size = int(args.target_size)

alive = set()
if args.symlink_roots_dir is not None:
    for path in Path(args.symlink_roots_dir).rglob("*"):
        if path.is_symlink():
            alive.add(path.resolve())

paths = []
total_size = 0
for cache_path in cache_dir.iterdir():
    stat = cache_path.stat()
    size = stat.st_size
    total_size += size
    age = time.time() - stat.st_atime
    score = size * age
    paths.append((cache_path, score, size))

print(f"Total cache size: {total_size / 1024. / 1024.:.1f} MiB")

paths.sort(key=lambda t: t[1])
while total_size >= target_size:
    cache_path, _, size = paths.pop()
    if cache_path in alive:
        print(f"{cache_path} is alive")
    else:
        cache_path.unlink(missing_ok=True)
        print(f"Removing {cache_path} of {size / 1024. / 1024.:.1f} MiB")
        total_size -= size
