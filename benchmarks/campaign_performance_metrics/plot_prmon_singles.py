#!/usr/bin/env python3
"""
Plot aggregate PRMON resource-usage histograms for all SINGLES datasets
in a given EPIC campaign.

Discovers datasets via:
    rucio list-dids --filter type=DATASET 'epic:/RECO/*<campaign>*SINGLE*'

For each discovered dataset a random sample of jobs is pulled from the
corresponding XRootD LOG directory, parsed, and pooled into a single
aggregate distribution.

Usage:
    python3 plot_prmon_singles.py --campaign 26.03.0 [--sample 5] [--seed 42]

    # Dry-run: just list the datasets that would be queried
    python3 plot_prmon_singles.py --campaign 26.03.0 --list-only
"""

import argparse
import subprocess
import random
import tarfile
import io
import csv
import re
from collections import defaultdict
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# ── Constants ────────────────────────────────────────────────────────────────

XROOT_BASE   = "root://dtn-eic.jlab.org"
RUCIO_SCOPE  = "epic"
RUCIO_PREFIX = "epic:/RECO"
LOG_PREFIX   = "/work/eic3/EPIC/LOGS/LOG"

STEPS = ["npsim", "eicrecon"]

PLOT_METRICS = [
    ("rss",   "Peak RSS Memory",  "MB", 1 / 1024),
    ("vmem",  "Peak Vmem",        "MB", 1 / 1024),
    ("wtime", "Wall Time",        "s",  1),
    ("utime", "CPU Time (user)",  "s",  1),
]


# ── Rucio helpers ─────────────────────────────────────────────────────────────

def rucio_list_singles(campaign):
    """
    Return a list of Rucio dataset DIDs matching the campaign SINGLES pattern.

    Runs:
        rucio list-dids --filter type=DATASET 'epic:/RECO/*<campaign>*SINGLE*'
    """
    pattern = f"{RUCIO_SCOPE}:/RECO/*{campaign}*SINGLE*"
    print(f"Running: rucio list-dids --filter type=DATASET '{pattern}'")
    result = subprocess.run(
        ["rucio", "list-dids", "--filter", "type=DATASET", pattern],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"[ERROR] rucio list-dids failed:\n{result.stderr.strip()}")
        return []

    datasets = []
    for line in result.stdout.splitlines():
        line = line.strip()
        # rucio output lines look like:  | epic:/RECO/... | DATASET |
        # or plain:  epic:/RECO/...
        m = re.search(r"(epic:/RECO/\S+)", line)
        if m:
            did = m.group(1).rstrip("|").strip()
            if "SINGLE" in did.upper() and campaign in did:
                datasets.append(did)

    return sorted(set(datasets))


# ── Path conversion ───────────────────────────────────────────────────────────

def rucio_to_log_dir(rucio_ds):
    """Convert a Rucio dataset name to an XRootD LOG directory path."""
    ds = rucio_ds.strip()
    if not ds.startswith(RUCIO_PREFIX):
        raise ValueError(
            f"Dataset does not start with '{RUCIO_PREFIX}': {ds!r}"
        )
    return LOG_PREFIX + ds[len(RUCIO_PREFIX):]


# ── XRootD helpers ────────────────────────────────────────────────────────────

def xrdfs_ls(remote_dir):
    result = subprocess.run(
        ["xrdfs", XROOT_BASE, "ls", remote_dir],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return []
    return [l.strip() for l in result.stdout.splitlines() if l.strip()]


def collect_tarballs(log_dir):
    """
    Walk up to two levels under log_dir to find *.log.tar.gz files.
    """
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
    result = subprocess.run(["xrdcp", "--silent", url, "-"], capture_output=True)
    return result.stdout if result.returncode == 0 else None


# ── Parsing ───────────────────────────────────────────────────────────────────

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
            summary[f"{k}_max"]  = max(vals)
            summary[f"{k}_mean"] = sum(vals) / len(vals)
    return summary


def parse_nevents_from_log(log_bytes):
    """Extract total event count from npsim or eicrecon log bytes."""
    text = log_bytes.decode("utf-8", errors="replace")
    m = re.search(r"Finished run \d+ after \d+ events \((\d+) events in total\)", text)
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
    """Returns {step: summary_dict}."""
    results = {}
    file_contents = {}
    try:
        with tarfile.open(fileobj=io.BytesIO(tarball_bytes), mode="r:gz") as tf:
            for member in tf.getmembers():
                f = tf.extractfile(member)
                if f:
                    file_contents[member.name] = f.read()
    except Exception as e:
        print(f"    [ERROR] {e}")
        return results

    for step in STEPS:
        prmon = next((c for n, c in file_contents.items() if f".{step}.prmon.txt" in n), None)
        if prmon is None:
            continue
        summary = parse_prmon_txt(prmon)
        if summary is None:
            continue
        log = next((c for n, c in file_contents.items() if f".{step}.log" in n), None)
        if log:
            n_ev = parse_nevents_from_log(log)
            if n_ev is not None:
                summary["n_events"] = n_ev
        results[step] = summary

    return results


# ── Data collection ───────────────────────────────────────────────────────────

def collect_data_for_dataset(rucio_ds, sample_size):
    """
    Returns {step: [per-job summary dicts]} for a single Rucio dataset.
    """
    log_dir = rucio_to_log_dir(rucio_ds)
    tarballs = collect_tarballs(log_dir)
    n = min(sample_size, len(tarballs))

    if not tarballs:
        print(f"    [WARN] No tarballs found in {log_dir}")
        return {step: [] for step in STEPS}

    print(f"    Found {len(tarballs)} tarballs, sampling {n}")
    sample = random.sample(tarballs, n)
    step_summaries = {step: [] for step in STEPS}

    for i, remote_path in enumerate(sample, 1):
        basename = remote_path.split("/")[-1]
        print(f"    [{i:2d}/{n}] {basename[:65]} ... ", end="", flush=True)
        data = stream_tarball_bytes(remote_path)
        if data is None:
            print("FAILED")
            continue
        results = extract_prmon_from_bytes(data)
        for step in STEPS:
            if step in results:
                step_summaries[step].append(results[step])
        found = [s for s in STEPS if s in results]
        print(f"ok ({', '.join(found)})")

    return step_summaries


def collect_all_singles(datasets, sample_size):
    """
    Iterate over all SINGLES datasets, collect samples, and return:
        aggregate: {step: [all per-job summaries]}
        per_ds:    {rucio_ds: {step: [per-job summaries]}}
    """
    aggregate  = {step: [] for step in STEPS}
    per_ds     = {}

    for idx, ds in enumerate(datasets, 1):
        short = ds.replace(RUCIO_PREFIX + "/", "")
        print(f"\n[{idx}/{len(datasets)}] {short}")
        ds_data = collect_data_for_dataset(ds, sample_size)
        per_ds[ds] = ds_data
        for step in STEPS:
            aggregate[step].extend(ds_data[step])

    return aggregate, per_ds


# ── Plotting ──────────────────────────────────────────────────────────────────

def make_aggregate_histograms(aggregate, campaign, n_datasets, timestamp):
    """
    aggregate : {step: [per-job summary dicts]}  (all datasets pooled)
    Produces one PNG per step with 2x2 metric histograms.
    """
    outfiles = []
    for step in STEPS:
        summaries = aggregate.get(step, [])
        n_jobs = len(summaries)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(
            f"PRMON aggregate — {step}  |  campaign: {campaign}\n"
            f"({n_jobs} jobs sampled across {n_datasets} SINGLES datasets)",
            fontsize=12, fontweight="bold"
        )
        axes = axes.flatten()

        for ax, (field, label, unit, scale) in zip(axes, PLOT_METRICS):
            key = f"{field}_max"
            vals = [s[key] * scale for s in summaries if key in s]

            if not vals:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", color="gray", fontsize=12)
                ax.set_title(label, fontsize=11)
                continue

            mu  = np.mean(vals)
            med = np.median(vals)
            p95 = np.percentile(vals, 95)

            bins = np.linspace(min(vals), max(vals), 30)
            ax.hist(vals, bins=bins, color="#4C72B0", alpha=0.75,
                    edgecolor="white", linewidth=0.5)
            ax.axvline(mu,  color="#DD8452", linestyle="--", linewidth=1.8,
                       label=f"mean {mu:.1f}")
            ax.axvline(med, color="#55A868", linestyle=":",  linewidth=1.8,
                       label=f"median {med:.1f}")
            ax.axvline(p95, color="#C44E52", linestyle="-.", linewidth=1.5,
                       label=f"p95 {p95:.1f}")

            ax.set_xlabel(f"{label} [{unit}]", fontsize=10)
            ax.set_ylabel("Jobs", fontsize=10)
            ax.set_title(label, fontsize=11)
            ax.legend(fontsize=8, framealpha=0.85)
            ax.grid(True, alpha=0.3, linestyle="--")

        plt.tight_layout()
        outfile = f"prmon_singles_{step}_{campaign}_{timestamp}.png"
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {outfile}")
        outfiles.append(outfile)

    return outfiles


def make_per_dataset_histograms(per_ds, campaign, timestamp):
    """
    One figure per step showing overlaid per-dataset distributions so outlier
    datasets are easy to spot. Each dataset gets its own color; legend shows
    the short dataset name and job count.
    """
    # Build a stable color cycle
    cmap = plt.get_cmap("tab20")
    ds_list = list(per_ds.keys())
    colors  = [cmap(i / max(len(ds_list), 1)) for i in range(len(ds_list))]

    outfiles = []
    for step in STEPS:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(
            f"PRMON per-dataset overlay — {step}  |  campaign: {campaign}",
            fontsize=12, fontweight="bold"
        )
        axes = axes.flatten()
        legend_handles = []

        for ax, (field, label, unit, scale) in zip(axes, PLOT_METRICS):
            key = f"{field}_max"

            # Compute global range for shared bins
            all_vals = []
            for ds in ds_list:
                summaries = per_ds[ds].get(step, [])
                all_vals.extend(s[key] * scale for s in summaries if key in s)

            if not all_vals:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", color="gray", fontsize=12)
                ax.set_title(label, fontsize=11)
                continue

            lo, hi = min(all_vals), max(all_vals)
            if lo == hi:
                lo -= 1; hi += 1
            bins = np.linspace(lo, hi, 20)

            for ci, ds in enumerate(ds_list):
                color     = colors[ci]
                summaries = per_ds[ds].get(step, [])
                vals      = [s[key] * scale for s in summaries if key in s]
                if not vals:
                    continue
                short = ds.replace(RUCIO_PREFIX + "/", "").split("/")[-1]
                ax.hist(vals, bins=bins, alpha=0.45, color=color,
                        edgecolor="white", linewidth=0.4)
                ax.axvline(np.mean(vals), color=color, linestyle="--",
                           linewidth=1.4, alpha=0.9)
                # Collect handles once (first metric subplot)
                if ax is axes[0]:
                    handle = mpatches.Patch(
                        color=color,
                        label=f"{short}  (n={len(vals)})"
                    )
                    legend_handles.append(handle)

            ax.set_xlabel(f"{label} [{unit}]", fontsize=10)
            ax.set_ylabel("Jobs", fontsize=10)
            ax.set_title(label, fontsize=11)
            ax.grid(True, alpha=0.3, linestyle="--")

        # Legend outside the subplots
        fig.legend(
            handles=legend_handles,
            loc="lower center",
            ncol=min(4, len(ds_list)),
            fontsize=7,
            framealpha=0.9,
            bbox_to_anchor=(0.5, 0.0),
        )
        n_rows = max(1, len(ds_list) // 4)
        plt.tight_layout(rect=[0, 0.04 * n_rows + 0.04, 1, 1])
        outfile = f"prmon_singles_perds_{step}_{campaign}_{timestamp}.png"
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {outfile}")
        outfiles.append(outfile)

    return outfiles


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate PRMON plots for all SINGLES datasets in a campaign."
    )
    parser.add_argument(
        "--campaign", required=True,
        help="Campaign string to match (e.g. '26.03.0').  "
             "Queries: rucio list-dids 'epic:/RECO/*<campaign>*SINGLE*'"
    )
    parser.add_argument(
        "--sample", type=int, default=5,
        help="Number of jobs to sample per dataset (default: 5)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--list-only", action="store_true",
        help="Just print the matching datasets and exit; do not fetch any data."
    )
    parser.add_argument(
        "--no-per-dataset", action="store_true",
        help="Skip the per-dataset overlay plot (only produce aggregate)."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    random.seed(args.seed)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # ── 1. Discover datasets ──────────────────────────────────────────────────
    print(f"\n=== Discovering SINGLES datasets for campaign: {args.campaign} ===")
    datasets = rucio_list_singles(args.campaign)

    if not datasets:
        print("[ERROR] No matching datasets found. Check campaign string and Rucio credentials.")
        raise SystemExit(1)

    print(f"\nFound {len(datasets)} dataset(s):")
    for ds in datasets:
        print(f"  {ds}")

    if args.list_only:
        raise SystemExit(0)

    # ── 2. Collect data ───────────────────────────────────────────────────────
    print(f"\n=== Collecting data ({args.sample} jobs/dataset) ===")
    aggregate, per_ds = collect_all_singles(datasets, args.sample)

    total_jobs = {step: len(aggregate[step]) for step in STEPS}
    print(f"\nTotal jobs collected: { {s: total_jobs[s] for s in STEPS} }")

    # ── 3. Plot ───────────────────────────────────────────────────────────────
    print("\n=== Generating plots ===")
    make_aggregate_histograms(aggregate, args.campaign, len(datasets), timestamp)

    if not args.no_per_dataset and len(datasets) > 1:
        make_per_dataset_histograms(per_ds, args.campaign, timestamp)

    print("\nDone.")
