#!/usr/bin/env python3
"""
Compare PRMON resource usage between two Rucio RECO datasets.

Rucio dataset names are converted to XRootD LOG paths by replacing the prefix:
    epic:/RECO  ->  /work/eic3/EPIC/LOGS/LOG

Usage:
    python3 plot_prmon.py \\
        --ds1 "epic:/RECO/26.03.0/epic_craterlake/Bkg_Exact1S_2us/GoldCt/10um/DIS/NC/10x100/minQ2=1" \\
        --ds2 "epic:/RECO/25.12.0/epic_craterlake/Bkg_Exactly1SignalPer2usFrame/GoldCoating/10um/DIS/NC/10x100/minQ2=1" \\
        [--label1 "26.03.0"] [--label2 "25.12.0"] [--sample 20]
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
RUCIO_PREFIX = "epic:/RECO"
LOG_PREFIX   = "/work/eic3/EPIC/LOGS/LOG"

STEPS = ["npsim", "eicrecon"]

PLOT_METRICS = [
    ("rss",   "Peak RSS Memory", "MB", 1/1024),
    ("vmem",  "Peak Vmem",       "MB", 1/1024),
    ("wtime", "Wall Time",       "s",  1),
    ("utime", "CPU Time (user)", "s",  1),
]

COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]


# ── Path conversion ──────────────────────────────────────────────────────────

def rucio_to_log_dir(rucio_ds):
    """Convert a Rucio dataset name to an XRootD LOG directory path."""
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
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return []
    return [l.strip() for l in result.stdout.splitlines() if l.strip()]


def collect_tarballs(log_dir):
    """
    The LOG dir may be structured as:
      log_dir/                     <- tarballs directly here, OR
      log_dir/minQ2=*/             <- one subdir level, OR
      log_dir/minQ2=*/             <- already included in log_dir
    Walk up to two levels to find *.log.tar.gz files.
    """
    tarballs = []
    top_entries = xrdfs_ls(log_dir)
    for entry in top_entries:
        if entry.endswith(".log.tar.gz"):
            tarballs.append(entry)          # flat layout
        else:
            for sub in xrdfs_ls(entry):     # one subdir level
                if sub.endswith(".log.tar.gz"):
                    tarballs.append(sub)
    return tarballs


def stream_tarball_bytes(remote_path):
    url = f"{XROOT_BASE}/{remote_path}"
    result = subprocess.run(["xrdcp", "--silent", url, "-"], capture_output=True)
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
            summary[f"{k}_max"]  = max(vals)
            summary[f"{k}_mean"] = sum(vals) / len(vals)
    return summary


def parse_nevents_from_log(log_bytes):
    """Extract total event count from npsim or eicrecon log bytes."""
    text = log_bytes.decode("utf-8", errors="replace")
    # npsim: 'Finished run 0 after N events (N events in total)'
    m = re.search(r"Finished run \d+ after \d+ events \((\d+) events in total\)", text)
    if m:
        return int(m.group(1))
    # eicrecon: 'Opened PODIO file "..." with N events'
    m = re.search(r'with (\d+) events \(format', text)
    if m:
        return int(m.group(1))
    # eicrecon fallback: last 'Status: N events processed'
    matches = re.findall(r"Status:\s+(\d+) events processed", text)
    if matches:
        return max(int(x) for x in matches)
    return None


def extract_prmon_from_bytes(tarball_bytes):
    """Returns {step: summary_dict}, summary includes 'n_events' if found in log."""
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


# ── Data collection ──────────────────────────────────────────────────────────

def collect_data(label, log_dir, sample_size):
    """Returns {step: [per-job summary dicts]}"""
    print(f"\n### {label}")
    print(f"    LOG dir: {log_dir}")
    tarballs = collect_tarballs(log_dir)
    n = min(sample_size, len(tarballs))
    print(f"    Found {len(tarballs)} tarballs, sampling {n}")

    if not tarballs:
        print("    [WARN] No tarballs found.")
        return {step: [] for step in STEPS}

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


# ── Plotting ─────────────────────────────────────────────────────────────────

def make_histograms(datasets, timestamp):
    """
    datasets: list of (label, rucio_ds, {step: [summaries]})
    One PNG per step, named with timestamp.
    A single figure-level legend identifies each dataset by its full Rucio ID.
    """
    outfiles = []
    for step in STEPS:
        fig, axes = plt.subplots(2, 2, figsize=(13, 9))
        fig.suptitle(f"PRMON distributions — {step}", fontsize=13, fontweight="bold")
        axes = axes.flatten()

        legend_handles = []

        for ax, (field, label, unit, scale) in zip(axes, PLOT_METRICS):
            key = f"{field}_max"

            all_vals = []
            for _, _ds, step_data in datasets:
                summaries = step_data.get(step, [])
                all_vals.extend(s[key] * scale for s in summaries if key in s)

            if not all_vals:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", color="gray", fontsize=12)
                ax.set_title(label, fontsize=11)
                continue

            lo, hi = min(all_vals), max(all_vals)
            if lo == hi:
                lo -= 1; hi += 1
            bins = np.linspace(lo, hi, 12)

            for ci, (clabel, rucio_ds, step_data) in enumerate(datasets):
                color = COLORS[ci % len(COLORS)]
                summaries = step_data.get(step, [])
                vals = [s[key] * scale for s in summaries if key in s]
                if not vals:
                    continue
                ev_counts = [s["n_events"] for s in summaries if "n_events" in s]
                ev_str = f", ~{int(np.median(ev_counts))} evts/job" if ev_counts else ""
                ax.hist(vals, bins=bins, alpha=0.60, color=color,
                        edgecolor="white", linewidth=0.5)
                ax.axvline(np.mean(vals), color=color, linestyle="--",
                           linewidth=1.8, alpha=0.95)
                # Build legend handle once (first subplot is enough)
                if len(legend_handles) < len(datasets):
                    handle = mpatches.Patch(
                        color=color,
                        label=f"{clabel}  (n={len(vals)}{ev_str})\n{rucio_ds}"
                    )
                    legend_handles.append(handle)

            ax.set_xlabel(f"{label} [{unit}]", fontsize=11)
            ax.set_ylabel("Jobs", fontsize=11)
            ax.set_title(label, fontsize=11)
            ax.grid(True, alpha=0.3, linestyle="--")

        # Single legend below all subplots
        fig.legend(
            handles=legend_handles,
            loc="lower center",
            ncol=len(datasets),
            fontsize=9,
            framealpha=0.9,
            bbox_to_anchor=(0.5, 0.0),
        )
        plt.tight_layout(rect=[0, 0.13 * len(datasets), 1, 1])
        outfile = f"prmon_compare_{step}_{timestamp}.png"
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {outfile}")
        outfiles.append(outfile)

    return outfiles


# ── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare PRMON logs for two Rucio RECO datasets."
    )
    parser.add_argument("--ds1", required=True,
                        help="First Rucio dataset name (epic:/RECO/...)")
    parser.add_argument("--ds2", required=True,
                        help="Second Rucio dataset name (epic:/RECO/...)")
    parser.add_argument("--label1", default=None,
                        help="Legend label for ds1 (default: auto from path)")
    parser.add_argument("--label2", default=None,
                        help="Legend label for ds2 (default: auto from path)")
    parser.add_argument("--sample", type=int, default=20,
                        help="Number of jobs to sample per dataset (default: 20)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    return parser.parse_args()


def auto_label(rucio_ds):
    """Generate a short label from the Rucio path (version + energy)."""
    parts = rucio_ds.replace("epic:/RECO/", "").split("/")
    # parts[0] = version, parts[-1] = energy or minQ2
    version = parts[0] if parts else "?"
    energy  = next((p for p in reversed(parts) if "x" in p and p[0].isdigit()), "")
    return f"{version} {energy}".strip()


if __name__ == "__main__":
    args = parse_args()
    random.seed(args.seed)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    log1 = rucio_to_log_dir(args.ds1)
    log2 = rucio_to_log_dir(args.ds2)

    label1 = args.label1 or auto_label(args.ds1)
    label2 = args.label2 or auto_label(args.ds2)

    print(f"\nDataset 1: {label1}")
    print(f"  Rucio : {args.ds1}")
    print(f"  LOG   : {log1}")
    print(f"\nDataset 2: {label2}")
    print(f"  Rucio : {args.ds2}")
    print(f"  LOG   : {log2}")

    data1 = collect_data(label1, log1, args.sample)
    data2 = collect_data(label2, log2, args.sample)

    print("\nGenerating histograms...")
    make_histograms([(label1, args.ds1, data1), (label2, args.ds2, data2)], timestamp)
    print("Done.")
