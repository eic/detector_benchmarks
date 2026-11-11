#!/usr/bin/env python3

import argparse
import os
import sys
import time
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Displays a text-based histogram of Snakemake cache entry size distribution in hourly bins",
    )
    parser.add_argument(
        "--cache-dir",
        "--snakemake-output-cache-dir",
        type=str,
        default=os.environ.get("SNAKEMAKE_OUTPUT_CACHE"),
        help="Path to Snakemake cache directory (defaults to $SNAKEMAKE_OUTPUT_CACHE)",
    )
    args = parser.parse_args()

    if not args.cache_dir or not os.path.isdir(args.cache_dir):
        return 0

    cache_dir = Path(args.cache_dir)
    try:
        entries = list(cache_dir.iterdir())
    except OSError:
        return 0

    now = time.time()
    bins = defaultdict(int)

    for p in entries:
        try:
            if not p.is_file():
                continue
            st = p.stat()
        except OSError:
            continue
        age_h = int(max(0, now - st.st_mtime) // 3600)
        bins[age_h] += st.st_size

    if not bins:
        print(f"Snakemake cache ({cache_dir}) is empty.")
        return 0

    max_h = max(bins)
    max_sz = max(bins.values())
    total_sz = sum(bins.values())
    print(f"Snakemake cache hourly size distribution ({total_sz / (1024**3):.1f} GiB total in {cache_dir}):")

    hours_to_print = range(max_h + 1) if max_h <= 72 else sorted(bins.keys())
    for h in hours_to_print:
        sz = bins[h]
        bar = "#" * int(40 * sz / max_sz) if max_sz else ""
        print(f"{h:3d}h: {sz / (1024**3):6.2f} GiB | {bar}".rstrip())

    return 0


if __name__ == "__main__":
    sys.exit(main())
