#!/usr/bin/env python3
"""
Aggregate per-dataset JSON files and produce PRMON histogram plots.

Reads all JSON files produced by collect_dataset.py, merges them, and
generates:
  1. Aggregate histograms (all datasets pooled)
  2. Per-dataset overlay histograms (one color per dataset)

Usage:
    python3 plot_singles_aggregate.py \
        --campaign 26.03.0 \
        --output-dir results/campaign_performance_metrics \
        collected/dataset_0.json collected/dataset_1.json ...
"""

import argparse
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

STEPS = ["npsim", "eicrecon"]

PLOT_METRICS = [
    ("rss",   "Peak RSS Memory",  "MB", 1 / 1024),
    ("vmem",  "Peak Vmem",        "MB", 1 / 1024),
    ("wtime", "Wall Time",        "s",  1),
    ("utime", "CPU Time (user)",  "s",  1),
]

RUCIO_PREFIX = "epic:/RECO"


# ── Data loading ─────────────────────────────────────────────────────────────

def load_and_merge(json_files):
    aggregate = {step: [] for step in STEPS}
    per_ds = {}

    for path in json_files:
        with open(path) as f:
            blob = json.load(f)
        ds_name = blob["dataset"]
        data = blob["data"]
        per_ds[ds_name] = data
        for step in STEPS:
            aggregate[step].extend(data.get(step, []))

    return aggregate, per_ds


# ── Plotting ─────────────────────────────────────────────────────────────────

def make_aggregate_histograms(aggregate, campaign, n_datasets, outdir):
    outfiles = []
    for step in STEPS:
        summaries = aggregate.get(step, [])
        n_jobs = len(summaries)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(
            f"PRMON aggregate \u2014 {step}  |  campaign: {campaign}\n"
            f"({n_jobs} jobs sampled across {n_datasets} SINGLES datasets)",
            fontsize=12, fontweight="bold",
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

            mu = np.mean(vals)
            med = np.median(vals)
            p95 = np.percentile(vals, 95)

            bins = np.linspace(min(vals), max(vals), 30)
            ax.hist(vals, bins=bins, color="#4C72B0", alpha=0.75,
                    edgecolor="white", linewidth=0.5)
            ax.axvline(mu, color="#DD8452", linestyle="--", linewidth=1.8,
                       label=f"mean {mu:.1f}")
            ax.axvline(med, color="#55A868", linestyle=":", linewidth=1.8,
                       label=f"median {med:.1f}")
            ax.axvline(p95, color="#C44E52", linestyle="-.", linewidth=1.5,
                       label=f"p95 {p95:.1f}")

            ax.set_xlabel(f"{label} [{unit}]", fontsize=10)
            ax.set_ylabel("Jobs", fontsize=10)
            ax.set_title(label, fontsize=11)
            ax.legend(fontsize=8, framealpha=0.85)
            ax.grid(True, alpha=0.3, linestyle="--")

        plt.tight_layout()
        outfile = f"{outdir}/prmon_singles_{step}_{campaign}.png"
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {outfile}")
        outfiles.append(outfile)

    return outfiles


def make_per_dataset_histograms(per_ds, campaign, outdir):
    import os

    ds_subdir = os.path.join(outdir, "per_dataset")
    os.makedirs(ds_subdir, exist_ok=True)

    step_colors = {"npsim": "#4C72B0", "eicrecon": "#DD8452"}

    outfiles = []
    for ds_name, data in per_ds.items():
        # Build a safe filename from the dataset path
        safe = ds_name.replace(RUCIO_PREFIX + "/", "").replace("/", "_")
        outfile = os.path.join(ds_subdir, f"{safe}.png")

        # 2 rows (steps) x 4 cols (metrics)
        fig, axes = plt.subplots(
            len(STEPS), len(PLOT_METRICS),
            figsize=(16, 5 * len(STEPS)),
            squeeze=False,
        )
        short_name = ds_name.replace(RUCIO_PREFIX + "/", "")
        fig.suptitle(
            f"PRMON \u2014 {short_name}\ncampaign: {campaign}",
            fontsize=11, fontweight="bold",
        )

        for row, step in enumerate(STEPS):
            summaries = data.get(step, [])
            for col, (field, label, unit, scale) in enumerate(PLOT_METRICS):
                ax = axes[row][col]
                key = f"{field}_max"
                vals = [s[key] * scale for s in summaries if key in s]

                ax.set_title(f"{step} \u2014 {label}", fontsize=9)
                ax.set_xlabel(f"{label} [{unit}]", fontsize=8)
                ax.set_ylabel("Jobs", fontsize=8)
                ax.grid(True, alpha=0.3, linestyle="--")

                if not vals:
                    ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                            ha="center", va="center", color="gray", fontsize=10)
                    continue

                mu  = np.mean(vals)
                med = np.median(vals)
                p95 = np.percentile(vals, 95)

                lo, hi = min(vals), max(vals)
                if lo == hi:
                    lo -= 0.5
                    hi += 0.5
                bins = np.linspace(lo, hi, max(5, min(20, len(vals))))

                color = step_colors.get(step, "#4C72B0")
                ax.hist(vals, bins=bins, color=color, alpha=0.75,
                        edgecolor="white", linewidth=0.5)
                ax.axvline(mu,  color="#55A868", linestyle="--", linewidth=1.6,
                           label=f"mean {mu:.1f}")
                ax.axvline(med, color="#C44E52", linestyle=":",  linewidth=1.6,
                           label=f"median {med:.1f}")
                ax.axvline(p95, color="#8172B2", linestyle="-.", linewidth=1.4,
                           label=f"p95 {p95:.1f}")
                ax.legend(fontsize=7, framealpha=0.85)

        plt.tight_layout()
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {outfile}")
        outfiles.append(outfile)

    return outfiles


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-dataset JSON files and produce plots.",
    )
    parser.add_argument("--campaign", required=True,
                        help="Campaign string (used in plot titles)")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for output plots")
    parser.add_argument("--no-per-dataset", action="store_true",
                        help="Skip per-dataset overlay plots")
    parser.add_argument("json_files", nargs="+",
                        help="Per-dataset JSON files from collect_dataset.py")
    args = parser.parse_args()

    print(f"Loading {len(args.json_files)} dataset file(s)...")
    aggregate, per_ds = load_and_merge(args.json_files)

    for step in STEPS:
        print(f"  {step}: {len(aggregate[step])} total jobs")

    print("\nGenerating aggregate plots...")
    make_aggregate_histograms(
        aggregate, args.campaign, len(per_ds), args.output_dir,
    )

    if not args.no_per_dataset and len(per_ds) > 1:
        print("Generating per-dataset overlay plots...")
        make_per_dataset_histograms(per_ds, args.campaign, args.output_dir)

    print("Done.")


if __name__ == "__main__":
    main()
