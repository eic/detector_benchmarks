"""EEEMCal reconstructed-cluster versus truth-cluster agreement.

The same event analysis supports electron particle guns and DIS. Particle-gun
events use MCParticles[0] and the highest-energy reconstructed cluster. DIS
events select the backward status-1 electron with the most negative pz from
MCScatteredElectrons and spatially match the nearest reconstructed cluster.

Residuals are filled only when the reference electron has exactly one
associated truth cluster. Truth-cluster multiplicity and energy fragmentation
are saved for every selected electron, so events with split truth clusters are
not silently discarded. These residuals measure reco--truth-cluster agreement,
not an independent intrinsic calorimeter resolution.
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


ENERGIES = (
    "100MeV", "200MeV", "500MeV", "1GeV",
    "2GeV", "5GeV", "10GeV", "20GeV",
)
RECO = "EcalEndcapNClusters"
TRUTH = "EcalEndcapNTruthClusters"
ASSOC = "_EcalEndcapNTruthClusterAssociations"
ETA_MIN = -3.5
ETA_MAX = -2.0


def read_file_list(path):
    paths = [
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not paths:
        raise ValueError(f"No input files in {path}")
    return paths


def reference_electron(event, sample_type):
    """Return (MCParticles index, eta, raw candidates, accepted candidates)."""
    if sample_type == "particle-gun":
        candidates = [0] if len(event["MCParticles.momentum.x"]) else []
    else:
        candidates = [
            int(index)
            for index in ak.to_list(event["MCScatteredElectrons_objIdx.index"])
        ]

    size = len(event["MCParticles.momentum.x"])
    if any(index < 0 or index >= size for index in candidates):
        raise ValueError("Scattered-electron relation points outside MCParticles")

    accepted = []
    for index in candidates:
        px, py, pz = (
            float(event[f"MCParticles.momentum.{axis}"][index]) for axis in "xyz"
        )
        pt = np.hypot(px, py)
        eta = np.arcsinh(pz / pt) if pt > 0 else np.nan
        if ETA_MIN < eta < ETA_MAX:
            accepted.append((pz, index, eta))
    if not accepted:
        return None, None, len(candidates), 0

    # This is the convention used by the existing EEEMCal DIS benchmark.
    _, index, eta = min(accepted)
    return index, eta, len(candidates), len(accepted)


def associated_truth_indices(event, particle_index):
    rec = ak.to_list(event[f"{ASSOC}_rec.index"])
    sim = ak.to_list(event[f"{ASSOC}_sim.index"])
    if len(rec) != len(sim):
        raise ValueError("Truth association rec/sim arrays differ in length")
    result = sorted({int(r) for r, s in zip(rec, sim) if int(s) == particle_index})
    if any(index < 0 or index >= len(event[f"{TRUTH}.energy"]) for index in result):
        raise ValueError("Truth association refers to an invalid cluster")
    return result


def select_reco_cluster(event, truth_index, sample_type, match_radius_mm):
    energies = np.asarray(ak.to_list(event[f"{RECO}.energy"]), dtype=float)
    if not len(energies):
        return None, None
    if sample_type == "particle-gun":
        return int(np.argmax(energies)), None

    tx = float(event[f"{TRUTH}.position.x"][truth_index])
    ty = float(event[f"{TRUTH}.position.y"][truth_index])
    rx = np.asarray(ak.to_list(event[f"{RECO}.position.x"]), dtype=float)
    ry = np.asarray(ak.to_list(event[f"{RECO}.position.y"]), dtype=float)
    distances = np.hypot(rx - tx, ry - ty)
    index = int(np.argmin(distances))
    distance = float(distances[index])
    return (index, distance) if distance < match_radius_mm else (None, distance)


def analyze(paths, sample_type, match_radius_mm):
    """Return strict residuals, fragmentation records, and selection counts."""
    counts = dict(
        events=0,
        no_reference_electron=0,
        multiple_reference_candidates=0,
        multiple_accepted_candidates=0,
        outside_pointing=0,
        zero_truth=0,
        multiple_truth=0,
        no_reco=0,
        reco_outside_match_radius=0,
        nonpositive_truth_energy=0,
        matched=0,
    )
    residuals = []
    fragmentation = []
    fields = ["MCParticles.momentum.*", f"{RECO}.*", f"{TRUTH}.*", f"{ASSOC}*"]
    if sample_type == "dis":
        fields.append("MCScatteredElectrons_objIdx.index")

    for events in uproot.iterate(
        {str(path): "events" for path in paths},
        filter_name=fields,
        step_size="100 MB",
    ):
        for event in events:
            counts["events"] += 1
            particle_index, eta, n_candidates, n_accepted = reference_electron(
                event, sample_type
            )
            if n_candidates > 1:
                counts["multiple_reference_candidates"] += 1
            if n_accepted > 1:
                counts["multiple_accepted_candidates"] += 1
            if not n_candidates:
                counts["no_reference_electron"] += 1
                continue
            if particle_index is None:
                counts["outside_pointing"] += 1
                continue

            truth_indices = associated_truth_indices(event, particle_index)
            truth_energies = np.asarray(
                [float(event[f"{TRUTH}.energy"][index]) for index in truth_indices])
            summed = float(np.sum(truth_energies)) if len(truth_energies) else 0.0
            largest = float(np.max(truth_energies)) if len(truth_energies) else 0.0
            fragmentation.append(
                (
                    eta,
                    len(truth_indices),
                    summed,
                    largest,
                    largest / summed if summed > 0 else np.nan,
                )
            )

            if len(truth_indices) != 1:
                counts["zero_truth" if not truth_indices else "multiple_truth"] += 1
                continue
            truth_index = truth_indices[0]
            truth_energy = float(event[f"{TRUTH}.energy"][truth_index])
            if truth_energy <= 0:
                counts["nonpositive_truth_energy"] += 1
                continue

            reco_index, match_distance = select_reco_cluster(
                event, truth_index, sample_type, match_radius_mm
            )
            if reco_index is None:
                key = "no_reco" if match_distance is None else "reco_outside_match_radius"
                counts[key] += 1
                continue
            dx, dy = (
                float(event[f"{RECO}.position.{axis}"][reco_index])
                - float(event[f"{TRUTH}.position.{axis}"][truth_index])
                for axis in "xy"
            )
            reco_energy = float(event[f"{RECO}.energy"][reco_index])
            residuals.append(
                ((reco_energy - truth_energy) / truth_energy, dx, dy, np.hypot(dx, dy))
            )
            counts["matched"] += 1

    return (
        np.asarray(residuals, dtype=float).reshape((-1, 4)),
        np.asarray(fragmentation, dtype=float).reshape((-1, 5)),
        counts,
    )


def central_width(values):
    if len(values) < 20:
        return None
    low, high = np.quantile(values, (0.16, 0.84))
    return float((high - low) / 2)


def plot_residuals(rows, label, output_path):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), layout="constrained")
    labels = (
        ("(E reco - E TC) / E TC", ""),
        ("Reco - TC x", "mm"),
        ("Reco - TC y", "mm"),
    )
    for i, (axis_label, unit) in enumerate(labels):
        axes[i].hist(rows[:, i], bins=80, histtype="step")
        axes[i].set_xlabel(f"{axis_label} [{unit}]" if unit else axis_label)
        axes[i].set_ylabel("Electron events")
    fig.suptitle(f"EEEMCal {label}: reco cluster versus unique associated TC")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_fragmentation(rows, label, output_path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.5), layout="constrained")
    multiplicity = rows[:, 1]
    max_multiplicity = int(np.max(multiplicity)) if len(rows) else 0
    axes[0].hist(
        multiplicity,
        bins=np.arange(-0.5, max(2.5, max_multiplicity + 1.5)),
        histtype="step",
    )
    axes[0].set_xlabel("Associated truth clusters")
    axes[0].set_ylabel("Selected electrons")
    axes[1].hist(
        rows[:, 4][np.isfinite(rows[:, 4])],
        bins=np.linspace(0, 1, 51),
        histtype="step",
    )
    axes[1].set_xlabel("Largest / summed TC energy")
    axes[1].set_ylabel("Selected electrons")
    axes[2].hist2d(
        rows[:, 0],
        rows[:, 1],
        bins=[35, np.arange(-0.5, max(2.5, max_multiplicity + 1.5))],
    )
    axes[2].set_xlabel("Generated electron eta")
    axes[2].set_ylabel("Associated truth clusters")
    fig.suptitle(f"EEEMCal {label}: truth-cluster fragmentation")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def result_summary(residuals, counts):
    return {
        "counts": counts,
        "central_68_half_width_energy_x_y": [
            central_width(residuals[:, i]) for i in range(3)
        ],
        "median_distance_mm": float(np.median(residuals[:, 3])) if len(residuals) else None,
    }


def run_particle_gun(args):
    summary = {}
    widths = []
    for energy in ENERGIES:
        paths = read_file_list(args.input_path_format.format(energy=energy))
        residuals, fragmentation, counts = analyze(
            paths, "particle-gun", args.match_radius_mm
        )
        np.savez_compressed(args.output_dir / f"residuals_{energy}.npz",
                            residuals=residuals, fragmentation=fragmentation)
        summary[energy] = result_summary(residuals, counts)
        widths.append(
            [
                np.nan if value is None else value
                for value in summary[energy]["central_68_half_width_energy_x_y"]
            ]
        )
        plot_residuals(residuals, energy, args.output_dir / f"residuals_{energy}.png")

    energy_values = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10, 20])
    widths = np.asarray(widths)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    axes[0].plot(energy_values, 100 * widths[:, 0], "o-")
    axes[0].set_ylabel("Energy residual central 68% half-width [%]")
    for index, label in ((1, "x"), (2, "y")):
        axes[1].plot(energy_values, widths[:, index], "o-", label=label)
    axes[1].set_ylabel("Position residual central 68% half-width [mm]")
    axes[1].legend()
    for axis in axes:
        axis.set_xlabel("Thrown electron energy [GeV]")
        axis.set_xscale("log")
    fig.savefig(args.output_dir / "agreement_widths.png", dpi=150)
    plt.close(fig)
    return summary


def run_dis(args):
    paths = read_file_list(args.input_file_list)
    residuals, fragmentation, counts = analyze(paths, "dis", args.match_radius_mm)
    np.savez_compressed(args.output_dir / "residuals.npz", residuals=residuals)
    np.savez_compressed(
        args.output_dir / "fragmentation.npz",
        eta=fragmentation[:, 0],
        multiplicity=fragmentation[:, 1],
        summed_energy=fragmentation[:, 2],
        largest_energy=fragmentation[:, 3],
        largest_fraction=fragmentation[:, 4],
    )
    plot_residuals(residuals, args.label, args.output_dir / "residuals.png")
    plot_fragmentation(fragmentation, args.label, args.output_dir / "fragmentation.png")
    return {"sample": result_summary(residuals, counts)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-type", choices=("particle-gun", "dis"), required=True)
    parser.add_argument("--input-path-format",
                        help="Particle-gun list-file pattern containing {energy}")
    parser.add_argument("--input-file-list", help="DIS reconstructed-file list")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--label", default="DIS")
    parser.add_argument("--match-radius-mm", type=float, default=50.0)
    args = parser.parse_args()
    if args.sample_type == "particle-gun" and not args.input_path_format:
        parser.error("--input-path-format is required for particle-gun inputs")
    if args.sample_type == "dis" and not args.input_file_list:
        parser.error("--input-file-list is required for DIS inputs")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = (
        run_particle_gun(args) if args.sample_type == "particle-gun" else run_dis(args)
    )
    with (args.output_dir / "summary.json").open("w") as stream:
        json.dump(
            {
                "sample_type": args.sample_type,
                "selection": (
                    "leading reco cluster; MCParticles[0]"
                    if args.sample_type == "particle-gun"
                    else "most-negative-pz accepted MCScatteredElectrons candidate; "
                    "nearest reco cluster"
                ),
                "strict_residual_selection": "exactly one associated truth cluster",
                "eta_range": [ETA_MIN, ETA_MAX],
                "match_radius_mm": (
                    args.match_radius_mm if args.sample_type == "dis" else None
                ),
                "interpretation": "reco-TC agreement, not intrinsic resolution",
                "results": results,
            },
            stream,
            indent=2,
        )


if __name__ == "__main__":
    main()
