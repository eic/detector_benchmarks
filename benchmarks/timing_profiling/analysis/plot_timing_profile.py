#!/usr/bin/env python3
"""
Analyze performance profile histograms to identify hot regions
Extracts summed timing values for regions in R and Z from the ZR plot
Creates annotated plots showing hot regions
"""

import sys
import os
import ROOT
import numpy as np

def analyze_hot_regions(input_file, output_prefix=None, n_events=None):
    """
    Analyze the m_zr histogram to find hot regions and create annotated plots

    Args:
        input_file: Path to the histos.root file
        output_prefix: Prefix for output files (default: input filename without .root)
        n_events: Number of events for normalization (default: 500)
    """
    # Set ROOT to batch mode
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)

    # Open the file
    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        print(f"Error: Cannot open file {input_file}")
        return 1

    # Get histograms
    m_zr = f.Get("m_zr")
    m_xy = f.Get("m_xy")

    if not m_zr:
        print("Error: Could not find histogram m_zr in file")
        f.Close()
        return 1

    if output_prefix is None:
        output_prefix = input_file.replace(".root", "")

    print(f"Analyzing hot regions from: {input_file}")
    print(f"=" * 80)
    print()

    # Number of events for normalization
    if n_events is None:
        n_events = 500

    print(f"Number of events (for normalization): {n_events}")

    # Get histogram properties
    nbins_z = m_zr.GetNbinsX()
    nbins_r = m_zr.GetNbinsY()
    z_min = m_zr.GetXaxis().GetXmin()
    z_max = m_zr.GetXaxis().GetXmax()
    r_min = m_zr.GetYaxis().GetXmin()
    r_max = m_zr.GetYaxis().GetXmax()

    print(f"Histogram dimensions:")
    print(f"  Z bins: {nbins_z}, range: [{z_min/1000:.2f}, {z_max/1000:.2f}] m")
    print(f"  R bins: {nbins_r}, range: [{r_min/1000:.2f}, {r_max/1000:.2f}] m")
    print()

    # Extract all bin values
    total_time = 0
    bin_values = []

    for iz in range(1, nbins_z + 1):
        for ir in range(1, nbins_r + 1):
            value = m_zr.GetBinContent(iz, ir)
            if value > 0:
                z_center = m_zr.GetXaxis().GetBinCenter(iz)
                r_center = m_zr.GetYaxis().GetBinCenter(ir)
                bin_values.append({
                    'z': z_center,
                    'r': r_center,
                    'time': value
                })
                total_time += value

    print(f"Total time: {total_time:.2e} ns = {total_time/1e9:.3f} s")
    print(f"Non-zero bins: {len(bin_values)} / {nbins_z * nbins_r}")
    print()

    if len(bin_values) == 0:
        print("No data in histogram!")
        f.Close()
        return 1

    # Define detector regions based on actual geometry
    # Format: (name, z_min, z_max, r_min, r_max)
    detector_regions = [
        ("Forward EM+Hadron Cal", 3346.25, 4971.75, 0, 2626),
        ("dRICH", 1944.75, 3250,  0, 1800),
        ("DIRC Bar (straight)", -2734, 1900, 750, 782.50),
        ("DIRC Bar (triangle)", -3150, -2734, 728.75, 1007.50),
        ("EEEMCal", -1945, -1740, 0, 540),
        ("Barrel Imaging and SciFi Cal",-2700, 1900,  826.77, 1217.61)
    ]

    # Group by region
    region_stats = {}

    for name, z_lo, z_hi, r_lo, r_hi in detector_regions:
        region_stats[name] = {
            'z_range': (z_lo, z_hi),
            'r_range': (r_lo, r_hi),
            'time': 0,
            'n_bins': 0
        }

    # Accumulate times by region
    for bin_data in bin_values:
        z = bin_data['z']
        r = bin_data['r']
        t = bin_data['time']

        for name, z_lo, z_hi, r_lo, r_hi in detector_regions:
            if z_lo <= z < z_hi and r_lo <= r < r_hi:
                region_stats[name]['time'] += t
                region_stats[name]['n_bins'] += 1
                break

    # Combine DIRC Bar straight and triangle sections for reporting
    dirc_combined_time = 0
    dirc_combined_bins = 0
    dirc_regions = []

    if "DIRC Bar (straight)" in region_stats:
        straight = region_stats["DIRC Bar (straight)"]
        dirc_combined_time += straight['time']
        dirc_combined_bins += straight['n_bins']
        dirc_regions.append(("DIRC Bar (straight)", straight))

    if "DIRC Bar (triangle)" in region_stats:
        triangle = region_stats["DIRC Bar (triangle)"]
        dirc_combined_time += triangle['time']
        dirc_combined_bins += triangle['n_bins']
        dirc_regions.append(("DIRC Bar (triangle)", triangle))

    # Create combined entry for sorting/display
    combined_stats = {}
    for name, stats in region_stats.items():
        if name not in ["DIRC Bar (straight)", "DIRC Bar (triangle)"]:
            combined_stats[name] = stats

    # Add combined DIRC entry
    if dirc_regions:
        combined_stats["DIRC"] = {
            'time': dirc_combined_time,
            'n_bins': dirc_combined_bins,
            'sub_regions': dirc_regions  # Keep track of sub-regions for drawing boxes
        }

    # Sort regions by time
    sorted_regions = sorted(combined_stats.items(),
                           key=lambda x: x[1]['time'],
                           reverse=True)

    # Print top regions
    print("Top regions by total time:")
    print("-" * 80)
    print(f"{'Region':<35} {'Time (ms)':<12} {'% Total':<10} {'Bins':<8}")
    print("-" * 80)

    for region_name, stats in sorted_regions:
        if stats['time'] > 0:
            time_ms = stats['time'] / 1e6
            percent = 100 * stats['time'] / total_time
            print(f"{region_name:<35} {time_ms:>10.2f}   {percent:>7.2f}%   {stats['n_bins']:>6}")

    print("-" * 80)
    print()

    # Create annotated ZR plot
    c_zr = ROOT.TCanvas("c_zr", "ZR Profile with Annotations", 1200, 800)
    c_zr.SetRightMargin(0.15)
    c_zr.SetLeftMargin(0.12)
    c_zr.SetLogz()

    m_zr.SetTitle(";Z [mm];R [mm];Time [ns]")
    m_zr.Draw("COLZ")

    # Draw region boundaries
    boxes = []
    labels = []

    # Draw boxes for each detector region
    for name, z_lo, z_hi, r_lo, r_hi in detector_regions:
        box = ROOT.TBox(z_lo, r_lo, z_hi, r_hi)
        box.SetFillStyle(0)
        box.SetLineColor(ROOT.kBlack)
        box.SetLineWidth(2)
        box.Draw()
        boxes.append(box)

    # Add text labels for top 5 hottest regions
    top_n = min(5, len([r for r in sorted_regions if r[1]['time'] > 0]))
    marker_boxes = []

    for idx, (region_name, stats) in enumerate(sorted_regions[:top_n]):
        if stats['time'] > 0:
            time_ms = stats['time'] / 1e6 / n_events  # Normalize by number of events, convert to ms
            percent = 100 * stats['time'] / total_time

            # Check if this is the combined DIRC
            if region_name == "DIRC" and 'sub_regions' in stats:
                # Draw thick boxes for both DIRC sub-regions
                all_z = []
                all_r = []
                for sub_name, sub_stats in stats['sub_regions']:
                    z_lo, z_hi = sub_stats['z_range']
                    r_lo, r_hi = sub_stats['r_range']

                    thick_box = ROOT.TBox(z_lo, r_lo, z_hi, r_hi)
                    thick_box.SetFillStyle(0)
                    thick_box.SetLineColor(ROOT.kBlack)
                    thick_box.SetLineWidth(4)
                    thick_box.Draw()
                    marker_boxes.append(thick_box)

                    all_z.extend([z_lo, z_hi])
                    all_r.extend([r_lo, r_hi])

                # Place label below the DIRC bar
                z_center = (min(all_z) + max(all_z)) / 2.0
                r_center = min(all_r) - 50  # Position below the region
                text_align = 23  # Bottom center alignment
            else:
                # Regular region
                z_lo, z_hi = stats['z_range']
                r_lo, r_hi = stats['r_range']

                # Calculate center of region
                z_center = (z_lo + z_hi) / 2.0
                r_center = (r_lo + r_hi) / 2.0

                # Draw thick black box outline for this region
                thick_box = ROOT.TBox(z_lo, r_lo, z_hi, r_hi)
                thick_box.SetFillStyle(0)
                thick_box.SetLineColor(ROOT.kBlack)
                thick_box.SetLineWidth(4)
                thick_box.Draw()
                marker_boxes.append(thick_box)

                text_align = 22  # Center alignment

            # Create label with integer values (time per event in ms)
            label = ROOT.TLatex(z_center, r_center, f"{idx+1}: {int(round(time_ms))}ms/evt ({int(round(percent))}%)")
            label.SetTextSize(0.025)
            label.SetTextAlign(text_align)
            label.SetTextColor(ROOT.kBlack)
            label.SetTextFont(62)  # Bold
            label.Draw()
            labels.append(label)

    # Add legend with top regions
    legend = ROOT.TLegend(0.12, 0.68, 0.50, 0.88)
    legend.SetTextSize(0.028)
    legend.SetTextColor(ROOT.kBlack)
    legend.SetFillStyle(1001)
    legend.SetFillColorAlpha(ROOT.kWhite, 0.85)
    legend.SetBorderSize(1)
    legend.SetMargin(0.05)  # Minimal left margin for entries
    legend.SetHeader("Most Computationally Expensive Regions", "C")

    for idx, (region_name, stats) in enumerate(sorted_regions[:top_n]):
        if stats['time'] > 0:
            time_ms = stats['time'] / 1e6 / n_events  # Normalize by number of events, convert to ms
            percent = 100 * stats['time'] / total_time
            legend.AddEntry(0, f"{idx+1}. {region_name}: {int(round(time_ms))}ms/evt ({int(round(percent))}%)", "")

    legend.Draw()

    c_zr.SaveAs(f"{output_prefix}_zr_annotated.png")
    c_zr.SaveAs(f"{output_prefix}_zr_annotated.pdf")
    print(f"Saved: {output_prefix}_zr_annotated.png")
    print(f"Saved: {output_prefix}_zr_annotated.pdf")

    # Create XY plot if available
    if m_xy:
        c_xy = ROOT.TCanvas("c_xy", "XY Profile", 800, 700)
        c_xy.SetRightMargin(0.15)
        c_xy.SetLogz()

        m_xy.SetTitle("Performance Profile - XY View;X [mm];Y [mm];Time [ns]")
        m_xy.Draw("COLZ")

        c_xy.SaveAs(f"{output_prefix}_xy.png")
        c_xy.SaveAs(f"{output_prefix}_xy.pdf")
        print(f"Saved: {output_prefix}_xy.png")
        print(f"Saved: {output_prefix}_xy.pdf")

    # Write CSV summary
    csv_file = f"{output_prefix}_regions.csv"
    with open(csv_file, 'w') as csvfile:
        csvfile.write("Region,Z_min(mm),Z_max(mm),R_min(mm),R_max(mm),Time(ns),Time(s),Percent,N_bins\n")
        for region_name, stats in sorted_regions:
            if stats['time'] > 0:
                time_ns = stats['time']
                time_s = time_ns / 1e9
                percent = 100 * time_ns / total_time
                n_bins = stats['n_bins']

                # Handle DIRC with sub-regions
                if region_name == "DIRC" and 'sub_regions' in stats:
                    # Calculate bounding box for all sub-regions
                    all_z = []
                    all_r = []
                    for sub_name, sub_stats in stats['sub_regions']:
                        z_lo, z_hi = sub_stats['z_range']
                        r_lo, r_hi = sub_stats['r_range']
                        all_z.extend([z_lo, z_hi])
                        all_r.extend([r_lo, r_hi])
                    z_lo, z_hi = min(all_z), max(all_z)
                    r_lo, r_hi = min(all_r), max(all_r)
                else:
                    z_lo, z_hi = stats['z_range']
                    r_lo, r_hi = stats['r_range']

                csvfile.write(f'"{region_name}",{z_lo},{z_hi},{r_lo},{r_hi},{time_ns:.2e},{time_s:.3f},{percent:.2f},{n_bins}\n')

    print(f"Saved: {csv_file}")
    print()

    f.Close()
    return 0

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: plot_timing_profile.py <histos.root> [output_prefix] [n_events]")
        print("\nExample:")
        print("  plot_timing_profile.py sim_output/timing_profile.histos.root results/timing_profiling/plot 500")
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else None
    n_events = int(sys.argv[3]) if len(sys.argv) > 3 else None

    sys.exit(analyze_hot_regions(input_file, output_prefix, n_events))
