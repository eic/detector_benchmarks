#!/usr/bin/env python3
"""Plot the histograms produced by the ``irt_analysis`` executable.

Several input files (e.g. one per generated particle species) can be given, in
which case the histograms are summed before plotting.
"""

import argparse
import os

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def load_histograms(input_files):
    """Sum the histograms of all input files, keyed by histogram name."""
    histograms = {}
    for input_file in input_files:
        f = ROOT.TFile.Open(input_file)
        if not f or f.IsZombie():
            raise OSError(f"Can not open {input_file}")
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            if not obj.InheritsFrom("TH1"):
                continue
            name = obj.GetName()
            if name in histograms:
                # TH1::Add checks binning consistency
                if not histograms[name].Add(obj):
                    raise ValueError(f"Inconsistent binning for histogram {name}")
            else:
                histograms[name] = obj.Clone(name)
        f.Close()
    return histograms


def save(canvas, path):
    canvas.SaveAs(path)
    canvas.Close()


def plot_th1(hist, path, title, xlabel):
    canvas = ROOT.TCanvas("c", "", 800, 600)
    hist.SetTitle(f"{title};{xlabel};Entries")
    hist.SetFillColorAlpha(ROOT.kAzure + 1, 0.7)
    hist.SetLineColor(ROOT.kAzure + 2)
    hist.SetStats(False)
    hist.Draw("HIST")
    save(canvas, path)


def plot_th2(hist, path, title, xlabel, ylabel):
    canvas = ROOT.TCanvas("c", "", 800, 600)
    canvas.SetRightMargin(0.15)
    if hist.GetMaximum() > 0:
        canvas.SetLogz()
    hist.SetTitle(f"{title};{xlabel};{ylabel};Entries")
    hist.SetStats(False)
    hist.Draw("COLZ")
    save(canvas, path)


def plot_confusion(hist, path, title):
    # normalize per true species (y axis) to get identification probabilities
    fractions = hist.Clone("Confusion_fractions")
    for j in range(1, hist.GetNbinsY() + 1):
        total = sum(hist.GetBinContent(i, j) for i in range(1, hist.GetNbinsX() + 1))
        for i in range(1, hist.GetNbinsX() + 1):
            fractions.SetBinContent(
                i, j, hist.GetBinContent(i, j) / total if total > 0 else 0.0
            )
            fractions.SetBinError(i, j, 0.0)

    canvas = ROOT.TCanvas("c", "", 800, 700)
    canvas.SetRightMargin(0.15)
    ROOT.gStyle.SetPaintTextFormat(".2f")
    fractions.SetTitle(f"{title};Reconstructed PID;True PID;Fraction of true species")
    fractions.SetStats(False)
    fractions.SetMinimum(0.0)
    fractions.SetMaximum(1.0)
    fractions.GetXaxis().SetLabelSize(0.06)
    fractions.GetYaxis().SetLabelSize(0.06)
    fractions.SetMarkerSize(1.5)
    fractions.Draw("COLZ TEXT")
    save(canvas, path)


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
            "p [GeV/c]",
        )
    if "m_RICH_ETA" in histograms:
        plot_th1(
            histograms["m_RICH_ETA"],
            out("eta.png"),
            title("Track pseudorapidity"),
            "#eta",
        )

    for prefix, filename, xlabel, text in [
        ("m_Theta_Rad_", "theta_rad_", "#theta_{C} [mrad]", "Cherenkov angle, radiator"),
        ("m_Npe_Rad_", "npe_rad_", "N_{pe}", "Number of photoelectrons, radiator"),
        ("m_NHits_Rad_", "nhits_rad_", "N_{hits}", "Number of hits, radiator"),
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
            "p [GeV/c]",
            "#theta_{C} [mrad]",
        )

    if "Confusion" in histograms:
        plot_confusion(
            histograms["Confusion"],
            out("confusion_matrix.png"),
            title("PID confusion matrix"),
        )


if __name__ == "__main__":
    main()
