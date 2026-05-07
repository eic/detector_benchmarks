#!/usr/bin/env python3
"""
Discover Rucio SINGLES datasets for a given EPIC campaign.

Writes one dataset DID per line to stdout (or a file via --output).

Usage:
    python3 discover_datasets.py --campaign 26.03.0
    python3 discover_datasets.py --campaign 26.03.0 --output datasets.txt
"""

import argparse
import re
import subprocess
import sys

RUCIO_SCOPE = "epic"
RUCIO_PREFIX = "epic:/RECO"


def rucio_list_singles(campaign):
    pattern = f"{RUCIO_SCOPE}:/RECO/*{campaign}*SINGLE*"
    print(f"Running: rucio list-dids --filter type=DATASET '{pattern}'",
          file=sys.stderr)
    result = subprocess.run(
        ["rucio", "list-dids", "--filter", "type=DATASET", pattern],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"[ERROR] rucio list-dids failed:\n{result.stderr.strip()}",
              file=sys.stderr)
        return []

    datasets = []
    for line in result.stdout.splitlines():
        line = line.strip()
        m = re.search(r"(epic:/RECO/\S+)", line)
        if m:
            did = m.group(1).rstrip("|").strip()
            if "SINGLE" in did.upper() and campaign in did:
                datasets.append(did)

    return sorted(set(datasets))


def main():
    parser = argparse.ArgumentParser(
        description="Discover SINGLES datasets for an EPIC campaign.",
    )
    parser.add_argument("--campaign", required=True,
                        help="Campaign string to match (e.g. '26.03.0')")
    parser.add_argument("--output", default=None,
                        help="Output file (default: stdout)")
    args = parser.parse_args()

    datasets = rucio_list_singles(args.campaign)

    if not datasets:
        print("[ERROR] No matching datasets found.", file=sys.stderr)
        raise SystemExit(1)

    print(f"Found {len(datasets)} dataset(s)", file=sys.stderr)

    if args.output:
        with open(args.output, "w") as f:
            for ds in datasets:
                f.write(ds + "\n")
        print(f"Written to {args.output}", file=sys.stderr)
    else:
        for ds in datasets:
            print(ds)


if __name__ == "__main__":
    main()
