#!/usr/bin/env python3
"""
Multi-Particle PFA Validation Analysis

Analyzes multi-particle event reconstruction performance for PandoraPFA, ArborPFA,
and baseline particle_flow algorithms.

Metrics:
- Particle multiplicity (reconstructed vs truth)
- Jet energy resolution (for di-jet events)
- Missing energy resolution
- Particle separation efficiency vs ΔR
- Event-level energy balance
"""

import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

try:
    import uproot
    import awkward as ak
except ImportError:
    print("ERROR: Required packages not found. Install with:")
    print("  pip install uproot awkward matplotlib")
    sys.exit(1)


def load_events(file_pattern):
    """Load events from ROOT files matching pattern"""
    files = glob.glob(file_pattern)
    if not files:
        raise FileNotFoundError(f"No files match pattern: {file_pattern}")
    
    print(f"Loading {len(files)} files...")
    events = []
    
    for fname in files:
        try:
            with uproot.open(fname) as f:
                tree = f["events"]
                events.append(tree.arrays())
        except Exception as e:
            print(f"WARNING: Failed to load {fname}: {e}")
    
    if not events:
        raise RuntimeError("No events loaded successfully")
    
    return ak.concatenate(events)


def plot_multiplicity(data, output_dir):
    """Plot reconstructed vs true particle multiplicity"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Particle Multiplicity Analysis', fontsize=14)
    
    ax = axes[0, 0]
    ax.set_title("Multiplicity Distribution")
    ax.set_xlabel("Number of Particles")
    ax.set_ylabel("Events")
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    ax = axes[0, 1]
    ax.set_title("Reco vs Truth Multiplicity")
    ax.set_xlabel("Truth Multiplicity")
    ax.set_ylabel("Reconstructed Multiplicity")
    ax.grid(True, alpha=0.3)
    # Diagonal line
    ax.plot([0, 20], [0, 20], 'k--', alpha=0.5, label='Ideal')
    
    ax = axes[1, 0]
    ax.set_title("Multiplicity Resolution")
    ax.set_xlabel("(N_reco - N_truth) / N_truth")
    ax.set_ylabel("Events")
    ax.grid(True, alpha=0.3)
    ax.axvline(0, color='black', linestyle='--', alpha=0.5)
    
    ax = axes[1, 1]
    ax.set_title("Algorithm Comparison")
    ax.set_xlabel("Algorithm")
    ax.set_ylabel("Multiplicity Resolution RMS [%]")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(output_dir) / "multiplicity_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_jet_resolution(data, output_dir):
    """Plot jet energy resolution for di-jet events"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Di-Jet Energy Resolution', fontsize=14)
    
    ax = axes[0]
    ax.set_title("Jet Energy Resolution")
    ax.set_xlabel("True Jet Energy [GeV]")
    ax.set_ylabel("(E_reco - E_true) / E_true")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    
    ax = axes[1]
    ax.set_title("Jet Mass Resolution")
    ax.set_xlabel("True Jet Mass [GeV]")
    ax.set_ylabel("(M_reco - M_true) / M_true")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    output_path = Path(output_dir) / "jet_energy_resolution.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_missing_energy(data, output_dir):
    """Plot missing energy resolution"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Missing Energy Analysis', fontsize=14)
    
    ax = axes[0]
    ax.set_title("Missing ET Resolution")
    ax.set_xlabel("True Missing ET [GeV]")
    ax.set_ylabel("Reconstructed Missing ET [GeV]")
    ax.grid(True, alpha=0.3)
    # Diagonal line
    ax.plot([0, 50], [0, 50], 'k--', alpha=0.5, label='Ideal')
    
    ax = axes[1]
    ax.set_title("Missing ET Distribution")
    ax.set_xlabel("(MET_reco - MET_true) [GeV]")
    ax.set_ylabel("Events")
    ax.grid(True, alpha=0.3)
    ax.axvline(0, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    output_path = Path(output_dir) / "missing_energy_resolution.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_energy_balance(data, output_dir):
    """Plot event-level energy balance"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_title('Event Energy Balance', fontsize=14)
    
    ax.set_xlabel("True Total Energy [GeV]")
    ax.set_ylabel("Reconstructed Total Energy [GeV]")
    ax.grid(True, alpha=0.3)
    # Diagonal line
    ax.plot([0, 100], [0, 100], 'k--', alpha=0.5, label='Ideal')
    ax.legend()
    
    plt.tight_layout()
    output_path = Path(output_dir) / "energy_balance.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_algorithm_performance(data, output_dir):
    """Create overall algorithm performance comparison"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_title('Multi-Particle Algorithm Performance', fontsize=14)
    
    # Placeholder data
    algorithms = ['PandoraPFA', 'ArborPFA', 'Baseline']
    metrics = {
        'Multiplicity\nResolution [%]': [15, 18, 22],
        'Energy\nBalance [%]': [8, 10, 12],
        'Missing ET\nResolution [%]': [20, 22, 25],
    }
    
    x = np.arange(len(algorithms))
    width = 0.25
    
    for i, (metric, values) in enumerate(metrics.items()):
        offset = (i - 1) * width
        ax.bar(x + offset, values, width, label=metric, alpha=0.8)
    
    ax.set_ylabel('Resolution / Error [%]')
    ax.set_xticks(x)
    ax.set_xticklabels(algorithms)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(output_dir) / "algorithm_performance.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="PFA Multi-Particle Analysis")
    parser.add_argument("--input-pattern", required=True,
                        help="Glob pattern for input ROOT files")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for plots")
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Analyzing multi-particle events...")
    print(f"Input pattern: {args.input_pattern}")
    print(f"Output directory: {args.output_dir}")
    
    try:
        # Load events
        events = load_events(args.input_pattern)
        print(f"Loaded {len(events)} events")
        
        # TODO: Extract and analyze collections
        data = {}  # Placeholder
        
        # Create plots
        print("\nGenerating plots...")
        plot_multiplicity(data, output_dir)
        plot_jet_resolution(data, output_dir)
        plot_missing_energy(data, output_dir)
        plot_energy_balance(data, output_dir)
        plot_algorithm_performance(data, output_dir)
        
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE")
        print("="*60)
        print(f"Results saved to: {output_dir}")
        print("\nNOTE: This is a placeholder implementation.")
        print("Production version requires actual data extraction and analysis.")
        
    except Exception as e:
        print(f"\nERROR: Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
