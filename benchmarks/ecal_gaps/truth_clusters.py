"""Truth-cluster integrity diagnostics; particle selection is caller-owned.

Association indices, not association row order, identify clusters and particles.
Multiple links to the same cluster are deduplicated. Cluster energy is counted
in full, not weighted by association weight: these are association diagnostics,
not a decomposition of shared-cluster energy into particle contributions.

The ecal_gaps workflow runs this analysis on reconstructed electron samples
and reports zero/one/multiple associated truth clusters versus generated eta.
Summed and largest associated-cluster energy responses and the largest-cluster
fraction are shown separately, so fragmentation is not hidden by an energy sum.
Numerical results are saved in JSON and per-sample NPZ files.

The reference particle is MCParticles[0] only for these particle-gun samples.
The association helper accepts an arbitrary MCParticles index, allowing future
DIS callers to supply the selected scattered electron. Missing associations are
retained as zero response. Outside a subsystem's acceptance, zero associations
are expected and should not be interpreted as inefficiency.

Backward/forward truth clusters and merged barrel truth clusters may use
different algorithms. The barrel output is labelled EcalBarrel, not ScFi or
imaging separately. These clean samples establish a baseline; they do not test
whether fragmentation under beam-background overlay has been fixed. No
truth-cluster position is used as a reference for angular resolution here.
"""

import argparse
import glob
import json
from pathlib import Path

import awkward as ak
import matplotlib
import numpy as np
import uproot

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def associated_cluster_energies(energies, rec_indices, sim_indices, particle_index):
    """Return distinct cluster energies associated with one MCParticles index.

    This function also accepts DIS scattered-electron indices selected by a
    caller. Invalid links fail explicitly rather than biasing the distributions.
    """
    if len(rec_indices) != len(sim_indices):
        raise ValueError("Association rec/sim arrays have different lengths")
    indices = sorted({r for r, s in zip(rec_indices, sim_indices) if s == particle_index})
    if any(r < 0 or r >= len(energies) for r in indices):
        raise ValueError("Association points outside the truth-cluster collection")
    return np.asarray([energies[r] for r in indices], dtype=float)


def analyze(paths, subsystem):
    collection = f"{subsystem}TruthClusters"
    association = f"_{subsystem}TruthClusterAssociations"
    branches = ["MCParticles.momentum.*", f"{collection}.energy",
                f"{association}_rec.index", f"{association}_sim.index"]
    records = []
    for events in uproot.iterate({p: "events" for p in paths},
                                 filter_name=branches, step_size="100 MB"):
        for event in events:
            # Only this selection is particle-gun-specific. DIS callers must
            # select a scattered electron instead of assuming MCParticles[0].
            particle_index = 0
            px, py, pz = [float(event[f"MCParticles.momentum.{c}"][particle_index])
                          for c in "xyz"]
            pt = np.hypot(px, py)
            momentum = np.hypot(pt, pz)
            if pt <= 0 or momentum <= 0:
                continue
            energies = associated_cluster_energies(
                ak.to_list(event[f"{collection}.energy"]),
                ak.to_list(event[f"{association}_rec.index"]),
                ak.to_list(event[f"{association}_sim.index"]), particle_index)
            total = float(np.sum(energies))
            largest = float(np.max(energies)) if len(energies) else 0.0
            records.append((np.arcsinh(pz / pt), len(energies), total / momentum,
                            largest / momentum, largest / total if total > 0 else np.nan))
    return np.asarray(records, dtype=float).reshape((-1, 5))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--detector-config", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    eta_edges = np.linspace(-4, 4, 41)
    centers = (eta_edges[:-1] + eta_edges[1:]) / 2
    results = {}
    for energy in ["500MeV", "5GeV", "20GeV"]:
        paths = sorted(glob.glob(
            f"sim_output/ecal_gaps/{args.detector_config}/e-/{energy}/*/*.eicrecon.edm4eic.root"))
        if not paths:
            raise RuntimeError(f"No reconstructed inputs for {energy}")
        for subsystem in ["EcalEndcapN", "EcalBarrel", "EcalEndcapP"]:
            data = analyze(paths, subsystem)
            key = f"{energy}_{subsystem}"
            np.savez_compressed(args.output_dir / f"truth_clusters_{key}.npz",
                                eta=data[:, 0], multiplicity=data[:, 1],
                                summed_response=data[:, 2], largest_response=data[:, 3],
                                largest_fraction=data[:, 4], eta_edges=eta_edges)
            counts = np.histogram(data[:, 0], eta_edges)[0]
            fractions = []
            means = []
            for low, high in zip(eta_edges[:-1], eta_edges[1:]):
                rows = data[(data[:, 0] >= low) & (data[:, 0] < high)]
                fractions.append([float(np.mean(condition)) if len(rows) else np.nan
                                  for condition in [rows[:, 1] == 0, rows[:, 1] == 1,
                                                    rows[:, 1] > 1]])
                means.append([float(np.mean(rows[:, c])) if len(rows) else np.nan
                              for c in [2, 3]])
            fractions, means = np.asarray(fractions), np.asarray(means)
            # Zero-association events are included; zero is not an efficiency
            # failure outside a subsystem's acceptance. Compare with hit response.
            results[key] = {"counts": counts.tolist(),
                            "association_fractions_zero_one_multiple":
                                [[float(v) if np.isfinite(v) else None for v in row]
                                 for row in fractions],
                            "mean_summed_and_largest_response":
                                [[float(v) if np.isfinite(v) else None for v in row]
                                 for row in means]}
            fig, axes = plt.subplots(1, 3, figsize=(13, 4))
            for i, label in enumerate(["zero", "one", "multiple"]):
                axes[0].plot(centers, fractions[:, i], label=label)
            axes[0].set_ylabel("Associated-cluster fraction")
            axes[0].set_ylim(0, 1.05)
            axes[0].legend()
            for i, label in enumerate(["sum", "largest"]):
                axes[1].plot(centers, means[:, i], label=label)
            axes[1].set_ylabel("Mean truth-cluster energy / thrown momentum")
            axes[1].legend()
            axes[2].hist2d(data[:, 0], data[:, 4], bins=[eta_edges, np.linspace(0, 1.01, 51)])
            axes[2].set_ylabel("Largest / summed associated energy")
            for axis in axes:
                axis.set_xlabel("Thrown electron eta")
            fig.suptitle(f"{energy}: {subsystem} truth-cluster diagnostics")
            fig.tight_layout()
            fig.savefig(args.output_dir / f"truth_clusters_{key}.png", dpi=150)
            plt.close(fig)
    with (args.output_dir / "truth_clusters_summary.json").open("w") as stream:
        json.dump({"eta_edges": eta_edges.tolist(), "results": results,
                   "energy_convention": "full energy of distinct associated clusters",
                   "particle_selection": "MCParticles index 0 (particle gun)"}, stream, indent=2)


if __name__ == "__main__":
    main()
