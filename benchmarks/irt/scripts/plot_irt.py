#!/usr/bin/env python3
"""Plot the histograms produced by the ``irt_analysis`` executable.

Several input files (e.g. one per generated particle species) can be given, in
which case the histograms are summed before plotting.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LogNorm

PID_LABELS = [r"$e$", r"$\pi$", r"$K$", r"$p$", r"$\mu$"]


def load_histograms(input_files):
    """Sum the histograms of all input files, keyed by histogram name."""
    histograms = {}
    for input_file in input_files:
        with uproot.open(input_file) as f:
            for key in f.keys():
                name = key.split(";")[0]
                values, *edges = f[key].to_numpy()
                if name in histograms:
                    previous_values, previous_edges = histograms[name]
                    if not all(
                        np.array_equal(a, b) for a, b in zip(previous_edges, edges)
                    ):
                        raise ValueError(f"Inconsistent binning for histogram {name}")
                    values = previous_values + values
                histograms[name] = (values, edges)
    return histograms


def plot_th1(hist, path, title, xlabel):
    values, (edges,) = hist
    fig, ax = plt.subplots()
    ax.stairs(values, edges, fill=True, alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Entries")
    ax.set_title(title)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_th2(hist, path, title, xlabel, ylabel):
    values, (xedges, yedges) = hist
    norm = LogNorm(vmin=1, vmax=max(values.max(), 1)) if values.max() > 0 else None
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(xedges, yedges, values.T, norm=norm)
    fig.colorbar(mesh, ax=ax, label="Entries")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_confusion(hist, path, title):
    values, _ = hist
    # normalize per true species (y axis) to get identification probabilities
    norm = values.sum(axis=0, keepdims=True)
    fractions = np.divide(values, norm, out=np.zeros_like(values), where=norm > 0)

    fig, ax = plt.subplots()
    mesh = ax.imshow(fractions.T, origin="lower", vmin=0.0, vmax=1.0, cmap="viridis")
    fig.colorbar(mesh, ax=ax, label="Fraction of true species")
    ax.set_xticks(range(len(PID_LABELS)), PID_LABELS)
    ax.set_yticks(range(len(PID_LABELS)), PID_LABELS)
    ax.set_xlabel("Reconstructed PID")
    ax.set_ylabel("True PID")
    ax.set_title(title)
    for i in range(fractions.shape[0]):
        for j in range(fractions.shape[1]):
            if norm[0, j] == 0:
                continue
            ax.text(
                i,
                j,
                f"{fractions[i, j]:.2f}",
                ha="center",
                va="center",
                color="white" if fractions[i, j] < 0.5 else "black",
                fontsize="small",
            )
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_files", nargs="+", help="histogram files to plot")
    parser.add_argument("-o", "--output-dir", required=True, help="output directory")
    parser.add_argument("-t", "--title", default="", help="prefix for plot titles")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    histograms = load_histograms(args.input_files)

    def out(name):
        return os.path.join(args.output_dir, name)

    def title(text):
        return f"{args.title}: {text}" if args.title else text

    if "m_RICH_Mom" in histograms:
        plot_th1(
            histograms["m_RICH_Mom"],
            out("momentum.png"),
            title("Track momentum"),
            r"$p$ [GeV/$c$]",
        )
    if "m_RICH_ETA" in histograms:
        plot_th1(
            histograms["m_RICH_ETA"],
            out("eta.png"),
            title("Track pseudorapidity"),
            r"$\eta$",
        )

    for prefix, filename, xlabel, text in [
        ("m_Theta_Rad_", "theta_rad_", r"$\theta_C$ [mrad]", "Cherenkov angle, radiator"),
        ("m_Npe_Rad_", "npe_rad_", r"$N_{pe}$", "Number of photoelectrons, radiator"),
        ("m_NHits_Rad_", "nhits_rad_", r"$N_{hits}$", "Number of hits, radiator"),
    ]:
        for name in sorted(k for k in histograms if k.startswith(prefix)):
            index = name[len(prefix):]
            plot_th1(
                histograms[name],
                out(f"{filename}{index}.png"),
                title(f"{text} {index}"),
                xlabel,
            )

    for name in sorted(k for k in histograms if k.startswith("m_ThetaMom_Rad_")):
        index = name[len("m_ThetaMom_Rad_"):]
        plot_th2(
            histograms[name],
            out(f"theta_vs_momentum_rad_{index}.png"),
            title(f"Cherenkov angle vs. momentum, radiator {index}"),
            r"$p$ [GeV/$c$]",
            r"$\theta_C$ [mrad]",
        )

    if "Confusion" in histograms:
        plot_confusion(
            histograms["Confusion"],
            out("confusion_matrix.png"),
            title("PID confusion matrix"),
        )


if __name__ == "__main__":
    main()
