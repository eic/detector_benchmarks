#!/usr/bin/env python3
"""Truth-matched K0S -> pi+ pi- performance of a secondary vertex collection.

A K0S counts as reconstructable if both of its pions are the best (highest-weight) MC match of
some reconstructed track, and as found if a two-track secondary vertex contains exactly those two
tracks. The output figure has three panels:
  1. vertex-finding efficiency (found / reconstructable) vs. true decay radius
  2. invariant mass of all two-track vertices with chi2/ndf < 3 and r > 2 mm (no truth used)
  3. invariant mass of the truth-matched (found) K0S vertices
"""

import argparse

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot

K0S_MASS = 0.4976  # GeV
RADIUS_BINS = np.array([0, 2, 5, 10, 20, 50, 100, 200, 500])  # mm
MASS_BINS = np.linspace(0.2, 1.0, 41)  # GeV


def best_mc_match(rec, sim, weight):
    """Map each reconstructed track index to its highest-weight MC particle index."""
    best = {}
    for r, s, w in zip(rec, sim, weight):
        if r not in best or w > best[r][1]:
            best[r] = (s, w)
    return {r: s for r, (s, w) in best.items()}


def process_file(path, vertices_collection, result):
    V = vertices_collection
    R = "ReconstructedChargedParticles"
    A = "ReconstructedChargedParticleAssociations"
    P = "MCParticles"
    branches = {
        "pdg": f"{P}/{P}.PDG",
        "end_x": f"{P}/{P}.endpoint.x",
        "end_y": f"{P}/{P}.endpoint.y",
        "dau_begin": f"{P}/{P}.daughters_begin",
        "dau_end": f"{P}/{P}.daughters_end",
        "dau_index": f"_{P}_daughters/_{P}_daughters.index",
        "px": f"{R}/{R}.momentum.x",
        "py": f"{R}/{R}.momentum.y",
        "pz": f"{R}/{R}.momentum.z",
        "energy": f"{R}/{R}.energy",
        "assoc_rec": f"_{A}_rec/_{A}_rec.index",
        "assoc_sim": f"_{A}_sim/_{A}_sim.index",
        "assoc_weight": f"{A}/{A}.weight",
        "vtx_x": f"{V}/{V}.position.x",
        "vtx_y": f"{V}/{V}.position.y",
        "vtx_chi2": f"{V}/{V}.chi2",
        "vtx_ndf": f"{V}/{V}.ndf",
        "vtx_begin": f"{V}/{V}.associatedParticles_begin",
        "vtx_end": f"{V}/{V}.associatedParticles_end",
        "vtx_index": f"_{V}_associatedParticles/_{V}_associatedParticles.index",
    }
    tree = uproot.open(path)["events"]
    arrays = tree.arrays(list(branches.values()), library="ak")
    events = {name: ak.to_list(arrays[branch]) for name, branch in branches.items()}

    for ev in range(len(events["pdg"])):
        e = {name: values[ev] for name, values in events.items()}
        rec_to_mc = best_mc_match(e["assoc_rec"], e["assoc_sim"], e["assoc_weight"])
        matched_mc = set(rec_to_mc.values())

        # Truth K0S -> pi+ pi-, keyed by the (unordered) pair of pion MC indices
        k0s = {}
        for i, pdg in enumerate(e["pdg"]):
            if pdg != 310:
                continue
            daughters = e["dau_index"][e["dau_begin"][i]:e["dau_end"][i]]
            if len(daughters) == 2 and all(abs(e["pdg"][j]) == 211 for j in daughters):
                r_true = np.hypot(e["end_x"][i], e["end_y"][i])
                if all(j in matched_mc for j in daughters):
                    k0s[frozenset(daughters)] = r_true
                    result["reconstructable_r"].append(r_true)

        px, py, pz, energy = e["px"], e["py"], e["pz"], e["energy"]
        found = set()
        for v in range(len(e["vtx_x"])):
            tracks = e["vtx_index"][e["vtx_begin"][v]:e["vtx_end"][v]]
            if len(tracks) != 2:
                continue
            # Vertices reference ReconstructedParticles, which MatchClusters fills with a copy of
            # ReconstructedChargedParticles first, in the same order, so the indices coincide
            t1, t2 = tracks
            if max(t1, t2) >= len(energy):
                continue
            mass2 = ((energy[t1] + energy[t2]) ** 2 - (px[t1] + px[t2]) ** 2
                     - (py[t1] + py[t2]) ** 2 - (pz[t1] + pz[t2]) ** 2)
            mass = np.sqrt(max(mass2, 0))
            r_vtx = np.hypot(e["vtx_x"][v], e["vtx_y"][v])
            chi2_ndf = e["vtx_chi2"][v] / (e["vtx_ndf"][v] + 1e-10)  # as in invariant_mass.py
            if chi2_ndf < 3 and r_vtx > 2:
                result["displaced_mass"].append(mass)

            pair = frozenset((rec_to_mc.get(t1), rec_to_mc.get(t2)))
            if pair in k0s and pair not in found:
                found.add(pair)
                result["found_r"].append(k0s[pair])
                result["found_mass"].append(mass)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", "-i", nargs="+", required=True, help="Input edm4eic ROOT file(s)")
    parser.add_argument("--output", "-o", required=True, help="Output plot filename")
    parser.add_argument("--vertices-collection", default="SecondaryVerticesAMVF",
                        help="Name of the vertices collection (default: SecondaryVerticesAMVF)")
    args = parser.parse_args()

    result = {"reconstructable_r": [], "found_r": [], "found_mass": [], "displaced_mass": []}
    for path in args.input:
        print(f"Loading {path}...")
        process_file(path, args.vertices_collection, result)
    result = {name: np.array(values) for name, values in result.items()}

    n_reco, _ = np.histogram(result["reconstructable_r"], RADIUS_BINS)
    n_found, _ = np.histogram(result["found_r"], RADIUS_BINS)
    print(f"\n=== K0S -> pi+ pi- found in {args.vertices_collection} ===")
    print(f"{'true r [mm]':>14} {'found':>6} / {'reco':<6} efficiency")
    for lo, hi, k, n in zip(RADIUS_BINS[:-1], RADIUS_BINS[1:], n_found, n_reco):
        eff = f"{k / n:.0%}" if n else "-"
        print(f"{f'[{lo}, {hi})':>14} {k:>6} / {n:<6} {eff}")
    print(f"{'total':>14} {len(result['found_r']):>6} / {len(result['reconstructable_r']):<6}")

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    has = n_reco > 0
    eff = np.divide(n_found, n_reco, out=np.zeros(len(n_reco)), where=has)
    err = np.sqrt(np.divide(eff * (1 - eff), n_reco, out=np.zeros(len(n_reco)), where=has))
    centers = 0.5 * (RADIUS_BINS[1:] + RADIUS_BINS[:-1])
    axs[0].errorbar(centers[has], eff[has], yerr=err[has], xerr=0.5 * np.diff(RADIUS_BINS)[has], fmt="o")
    axs[0].set(xscale="symlog", xlim=(0, RADIUS_BINS[-1]), ylim=(0, 1.05),
               xlabel="True K⁰_S decay radius [mm]", ylabel="Vertex-finding efficiency",
               title=f"K⁰_S (both pions reconstructed) found in {args.vertices_collection}")

    axs[1].hist(result["displaced_mass"], bins=MASS_BINS, histtype="step", linewidth=2)
    axs[1].set(xlabel="Invariant mass [GeV]", ylabel="Vertices",
               title=f"{args.vertices_collection}: χ²/ndf < 3, r > 2 mm")

    axs[2].hist(result["found_mass"], bins=MASS_BINS, histtype="step", linewidth=2)
    axs[2].set(xlabel="Invariant mass [GeV]", ylabel="Vertices", title="Truth-matched K⁰_S vertices")

    for ax in axs[1:]:
        ax.axvline(K0S_MASS, color="red", linestyle="--", label="K⁰_S")
        ax.legend()

    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
