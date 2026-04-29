#!/usr/bin/env python3
"""
Collect PRMON data for a single Rucio dataset and write results as JSON.

Designed to be called by Snakemake as an independent job per dataset,
enabling parallel data collection across datasets.

Usage:
    python3 collect_dataset.py \
        --dataset "epic:/RECO/26.03.0/.../SINGLE/..." \
        --sample 5 --seed 42 \
        --output collected/dataset_0.json
"""

import argparse
import csv
import io
import json
import random
import re
import subprocess
import sys
import tarfile
from collections import defaultdict

XROOT_BASE = "root://dtn-eic.jlab.org"
RUCIO_PREFIX = "epic:/RECO"
LOG_PREFIX = "/work/eic3/EPIC/LOGS/LOG"

STEPS = ["npsim", "eicrecon"]


# ── Path conversion ──────────────────────────────────────────────────────────

def rucio_to_log_dir(rucio_ds):
    ds = rucio_ds.strip()
    if not ds.startswith(RUCIO_PREFIX):
        raise ValueError(
            f"Dataset does not start with '{RUCIO_PREFIX}': {ds!r}"
        )
    return LOG_PREFIX + ds[len(RUCIO_PREFIX):]


# ── XRootD helpers ───────────────────────────────────────────────────────────

def xrdfs_ls(remote_dir):
    result = subprocess.run(
        ["xrdfs", XROOT_BASE, "ls", remote_dir],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return []
    return [l.strip() for l in result.stdout.splitlines() if l.strip()]


def collect_tarballs(log_dir):
    tarballs = []
    top_entries = xrdfs_ls(log_dir)
    for entry in top_entries:
        if entry.endswith(".log.tar.gz"):
            tarballs.append(entry)
        else:
            for sub in xrdfs_ls(entry):
                if sub.endswith(".log.tar.gz"):
                    tarballs.append(sub)
    return tarballs


def stream_tarball_bytes(remote_path):
    url = f"{XROOT_BASE}/{remote_path}"
    result = subprocess.run(["xrdcp", "--silent", url, "-"],
                            capture_output=True)
    return result.stdout if result.returncode == 0 else None


# ── Parsing ──────────────────────────────────────────────────────────────────

def parse_prmon_txt(content_bytes):
    text = content_bytes.decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    rows = list(reader)
    if not rows:
        return None
    columns = defaultdict(list)
    for row in rows:
        for k, v in row.items():
            if k and k != "Time":
                try:
                    columns[k].append(float(v))
                except (ValueError, TypeError):
                    pass
    summary = {"n_records": len(rows)}
    for k, vals in columns.items():
        if vals:
            summary[f"{k}_max"] = max(vals)
            summary[f"{k}_mean"] = sum(vals) / len(vals)
    return summary


def parse_nevents_from_log(log_bytes):
    text = log_bytes.decode("utf-8", errors="replace")
    m = re.search(
        r"Finished run \d+ after \d+ events \((\d+) events in total\)", text
    )
    if m:
        return int(m.group(1))
    m = re.search(r'with (\d+) events \(format', text)
    if m:
        return int(m.group(1))
    matches = re.findall(r"Status:\s+(\d+) events processed", text)
    if matches:
        return max(int(x) for x in matches)
    return None


def extract_prmon_from_bytes(tarball_bytes):
    results = {}
    file_contents = {}
    try:
        with tarfile.open(fileobj=io.BytesIO(tarball_bytes),
                          mode="r:gz") as tf:
            for member in tf.getmembers():
                f = tf.extractfile(member)
                if f:
                    file_contents[member.name] = f.read()
    except Exception as e:
        print(f"    [ERROR] {e}", file=sys.stderr)
        return results

    for step in STEPS:
        prmon = next(
            (c for n, c in file_contents.items()
             if f".{step}.prmon.txt" in n),
            None,
        )
        if prmon is None:
            continue
        summary = parse_prmon_txt(prmon)
        if summary is None:
            continue
        log = next(
            (c for n, c in file_contents.items() if f".{step}.log" in n),
            None,
        )
        if log:
            n_ev = parse_nevents_from_log(log)
            if n_ev is not None:
                summary["n_events"] = n_ev
        results[step] = summary

    return results


# ── Data collection ──────────────────────────────────────────────────────────

def collect_data_for_dataset(rucio_ds, sample_size):
    log_dir = rucio_to_log_dir(rucio_ds)
    tarballs = collect_tarballs(log_dir)
    n = min(sample_size, len(tarballs))

    if not tarballs:
        print(f"    [WARN] No tarballs found in {log_dir}", file=sys.stderr)
        return {step: [] for step in STEPS}

    print(f"    Found {len(tarballs)} tarballs, sampling {n}",
          file=sys.stderr)
    sample = random.sample(tarballs, n)
    step_summaries = {step: [] for step in STEPS}

    for i, remote_path in enumerate(sample, 1):
        basename = remote_path.split("/")[-1]
        print(f"    [{i:2d}/{n}] {basename[:65]} ... ",
              end="", flush=True, file=sys.stderr)
        data = stream_tarball_bytes(remote_path)
        if data is None:
            print("FAILED", file=sys.stderr)
            continue
        results = extract_prmon_from_bytes(data)
        for step in STEPS:
            if step in results:
                step_summaries[step].append(results[step])
        found = [s for s in STEPS if s in results]
        print(f"ok ({', '.join(found)})", file=sys.stderr)

    return step_summaries


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Collect PRMON data for a single Rucio dataset.",
    )
    parser.add_argument("--dataset", required=True,
                        help="Rucio dataset DID (epic:/RECO/...)")
    parser.add_argument("--sample", type=int, default=5,
                        help="Number of jobs to sample (default: 5)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--output", required=True,
                        help="Output JSON file path")
    args = parser.parse_args()

    random.seed(args.seed)

    short = args.dataset.replace(RUCIO_PREFIX + "/", "")
    print(f"Collecting data for: {short}", file=sys.stderr)

    step_summaries = collect_data_for_dataset(args.dataset, args.sample)

    result = {
        "dataset": args.dataset,
        "sample_size": args.sample,
        "seed": args.seed,
        "data": step_summaries,
    }

    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)

    for step in STEPS:
        n = len(step_summaries[step])
        print(f"    {step}: {n} jobs collected", file=sys.stderr)

    print(f"Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
