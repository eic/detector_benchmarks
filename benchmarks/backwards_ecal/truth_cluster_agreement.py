"""EEEMCal reconstructed-cluster versus truth-cluster agreement.

This is separate from the endpoint-based position resolution in
backwards_ecal.org. It uses that benchmark's highest-energy reconstructed
EcalEndcapN cluster, requires one truth cluster associated with MCParticles[0],
and makes no spatial matching cut. The two clusters share reconstructed hits,
so their residuals are not an independent intrinsic position resolution.

Only electron particle-gun inputs are used. Counts of rejected events and
uncut residual arrays accompany the plots so long tails remain visible.
"""

import argparse
import json
from pathlib import Path

import awkward as ak
import matplotlib
import numpy as np
import uproot

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ENERGIES = ("100MeV", "200MeV", "500MeV", "1GeV", "2GeV", "5GeV", "10GeV", "20GeV")
RECO = "EcalEndcapNClusters"
TRUTH = "EcalEndcapNTruthClusters"
ASSOC = "_EcalEndcapNTruthClusterAssociations"


def associated_truth_indices(event):
    rec = ak.to_list(event[f"{ASSOC}_rec.index"])
    sim = ak.to_list(event[f"{ASSOC}_sim.index"])
    if len(rec) != len(sim):
        raise ValueError("Truth association rec/sim arrays differ in length")
    result = sorted({r for r, s in zip(rec, sim) if s == 0})
    if any(r < 0 or r >= len(event[f"{TRUTH}.energy"]) for r in result):
        raise ValueError("Truth association refers to an invalid cluster")
    return result


def analyze(paths):
    """Return uncut (energy fraction, dx, dy, distance) and selection counts."""
    counts = dict(events=0, outside_pointing=0, zero_truth=0, multiple_truth=0,
                  no_reco=0, nonpositive_truth_energy=0, matched=0)
    rows = []
    fields = ["MCParticles.momentum.*", f"{RECO}.*", f"{TRUTH}.*", f"{ASSOC}*"]
    for events in uproot.iterate({str(path): "events" for path in paths},
                                 filter_name=fields, step_size="100 MB"):
        for event in events:
            counts["events"] += 1
            px, py, pz = (float(event[f"MCParticles.momentum.{c}"][0]) for c in "xyz")
            pt = np.hypot(px, py)
            eta = np.arcsinh(pz / pt) if pt > 0 else np.nan
            if not -3.5 < eta < -2.0:
                counts["outside_pointing"] += 1
                continue
            truth_indices = associated_truth_indices(event)
            if len(truth_indices) != 1:
                counts["zero_truth" if not truth_indices else "multiple_truth"] += 1
                continue
            reco_energy = ak.to_list(event[f"{RECO}.energy"])
            if not reco_energy:
                counts["no_reco"] += 1
                continue
            t = truth_indices[0]
            et = float(event[f"{TRUTH}.energy"][t])
            if et <= 0:
                counts["nonpositive_truth_energy"] += 1
                continue
            r = int(np.argmax(reco_energy))
            dx, dy = (float(event[f"{RECO}.position.{c}"][r]) -
                      float(event[f"{TRUTH}.position.{c}"][t]) for c in "xy")
            rows.append(((float(reco_energy[r]) - et) / et, dx, dy, np.hypot(dx, dy)))
            counts["matched"] += 1
    return np.asarray(rows, dtype=float).reshape((-1, 4)), counts


def central_width(values):
    if len(values) < 20:
        return None
    low, high = np.quantile(values, (0.16, 0.84))
    return float((high - low) / 2)


def plot_residuals(rows, energy, output_dir):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), layout="constrained")
    for i, (label, unit) in enumerate((("(E reco - E TC) / E TC", ""),
                                       ("Reco - TC x", "mm"),
                                       ("Reco - TC y", "mm"))):
        axes[i].hist(rows[:, i], bins=80, histtype="step")
        axes[i].set_xlabel(f"{label} [{unit}]" if unit else label)
        axes[i].set_ylabel("Electron events")
    fig.suptitle(f"EEEMCal {energy}: leading reco cluster versus associated TC")
    fig.savefig(output_dir / f"residuals_{energy}.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-path-format", required=True,
                        help="List-file pattern containing {energy}; electron inputs only")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = {}
    widths = []
    for energy in ENERGIES:
        list_path = Path(args.input_path_format.format(energy=energy))
        paths = [Path(line.strip()) for line in list_path.read_text().splitlines() if line.strip()]
        if not paths:
            raise ValueError(f"No input files in {list_path}")
        rows, counts = analyze(paths)
        np.savez_compressed(args.output_dir / f"residuals_{energy}.npz", residuals=rows)
        width = [central_width(rows[:, i]) for i in range(3)]
        summary[energy] = {"counts": counts, "central_68_half_width_energy_x_y": width,
                           "median_distance_mm": float(np.median(rows[:, 3])) if len(rows) else None}
        widths.append([np.nan if v is None else v for v in width])
        plot_residuals(rows, energy, args.output_dir)
    energy_values = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10, 20])
    widths = np.asarray(widths)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    axes[0].plot(energy_values, 100 * widths[:, 0], "o-")
    axes[0].set_ylabel("Energy residual central 68% half-width [%]")
    for i, label in ((1, "x"), (2, "y")):
        axes[1].plot(energy_values, widths[:, i], "o-", label=label)
    axes[1].set_ylabel("Position residual central 68% half-width [mm]")
    axes[1].legend()
    for axis in axes:
        axis.set_xlabel("Thrown electron energy [GeV]")
        axis.set_xscale("log")
    fig.savefig(args.output_dir / "agreement_widths.png", dpi=150)
    plt.close(fig)
    with (args.output_dir / "summary.json").open("w") as stream:
        json.dump({"selection": "leading reco cluster; unique TC associated with MCParticles[0]; -3.5 < generated eta < -2.0; no distance cut",
                   "interpretation": "reco-TC agreement, not intrinsic resolution",
                   "results": summary}, stream, indent=2)


if __name__ == "__main__":
    main()
