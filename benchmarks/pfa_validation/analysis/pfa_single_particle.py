#!/usr/bin/env python3
"""
Single Particle PFA Validation Analysis

Analyzes single particle reconstruction performance for PandoraPFA, ArborPFA,
and baseline particle_flow algorithms.

Metrics:
- Energy resolution vs energy
- Momentum resolution vs momentum
- Angular resolution (theta, phi)
- Reconstruction efficiency vs energy and eta
- Particle ID accuracy
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


def extract_collections(events):
    """
    Extract particle collections from events
    
    TODO: Update collection names once confirmed by EICrecon session
    Expected collections:
    - MCParticles (truth)
    - PandoraPFAParticles (Pandora reconstruction)
    - ArborPFAParticles (Arbor reconstruction)
    - ReconstructedParticles or similar (baseline)
    """
    data = {}
    
    # Truth particles
    if "MCParticles" in events.fields:
        data["truth"] = events.MCParticles
    else:
        raise KeyError("MCParticles collection not found in events")
    
    # PandoraPFA (TODO: confirm exact collection name)
    pfa_names = [
        "PandoraPFAParticles",
        "PandoraPFAReconstructedParticles",
        "PandoraParticles"
    ]
    for name in pfa_names:
        if name in events.fields:
            data["pandora"] = events[name]
            break
    
    # ArborPFA (TODO: confirm exact collection name)
    arbor_names = [
        "ArborPFAParticles",
        "ArborPFAReconstructedParticles",
        "ArborParticles"
    ]
    for name in arbor_names:
        if name in events.fields:
            data["arbor"] = events[name]
            break
    
    # Baseline (TODO: confirm exact collection name)
    baseline_names = [
        "ReconstructedParticles",
        "ReconstructedChargedParticles",
        "ParticleFlowParticles"
    ]
    for name in baseline_names:
        if name in events.fields:
            data["baseline"] = events[name]
            break
    
    print(f"Found collections: {list(data.keys())}")
    return data


def match_truth_to_reco(truth, reco, dr_max=0.1):
    """
    Match truth particles to reconstructed particles by proximity
    
    Returns indices of matched reconstructed particles (-1 if no match)
    """
    # TODO: Implement proper truth matching
    # For now, use simple Delta R matching
    # Better approach: use truth links if available in EDM4eic
    
    # Placeholder implementation
    n_truth = len(truth)
    matched_indices = np.full(n_truth, -1, dtype=int)
    
    # This is a simplified matching - production code should use truth links
    print(f"WARNING: Using simplified truth matching. Production code should use EDM4eic truth links.")
    
    return matched_indices


def compute_energy_resolution(truth_e, reco_e):
    """Compute energy resolution (E_reco - E_truth) / E_truth"""
    return (reco_e - truth_e) / truth_e


def compute_efficiency(truth, reco, energy_bins):
    """Compute reconstruction efficiency vs energy"""
    # TODO: Implement efficiency calculation
    # Efficiency = N_reconstructed / N_truth in each energy bin
    pass


def plot_energy_resolution(data, particle_name, output_dir):
    """
    Plot energy resolution for all algorithms
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Energy Resolution - {particle_name}', fontsize=14)
    
    algorithms = ["pandora", "arbor", "baseline"]
    colors = {"pandora": "blue", "arbor": "red", "baseline": "green"}
    
    # TODO: Extract actual data and create plots
    # For now, create placeholder plots
    
    ax = axes[0, 0]
    ax.set_title("Energy Resolution vs True Energy")
    ax.set_xlabel("True Energy [GeV]")
    ax.set_ylabel("(E_reco - E_true) / E_true")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    
    ax = axes[0, 1]
    ax.set_title("Energy Resolution Distribution")
    ax.set_xlabel("(E_reco - E_true) / E_true")
    ax.set_ylabel("Counts")
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 0]
    ax.set_title("Reconstructed vs True Energy")
    ax.set_xlabel("True Energy [GeV]")
    ax.set_ylabel("Reconstructed Energy [GeV]")
    ax.grid(True, alpha=0.3)
    # Diagonal line
    ax.plot([0, 100], [0, 100], 'k--', alpha=0.5, label='Ideal')
    
    ax = axes[1, 1]
    ax.set_title("Algorithm Comparison")
    ax.set_xlabel("True Energy [GeV]")
    ax.set_ylabel("Energy Resolution RMS [%]")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(output_dir) / f"energy_resolution_vs_energy.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_efficiency(data, particle_name, output_dir):
    """Plot reconstruction efficiency"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f'Reconstruction Efficiency - {particle_name}', fontsize=14)
    
    ax = axes[0]
    ax.set_title("Efficiency vs Energy")
    ax.set_xlabel("True Energy [GeV]")
    ax.set_ylabel("Efficiency")
    ax.set_ylim([0, 1.1])
    ax.grid(True, alpha=0.3)
    ax.axhline(0.9, color='gray', linestyle='--', alpha=0.5, label='90% threshold')
    
    ax = axes[1]
    ax.set_title("Efficiency vs η")
    ax.set_xlabel("Pseudorapidity η")
    ax.set_ylabel("Efficiency")
    ax.set_ylim([0, 1.1])
    ax.grid(True, alpha=0.3)
    ax.axhline(0.9, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    output_path = Path(output_dir) / f"efficiency_vs_energy.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_algorithm_comparison(data, particle_name, output_dir):
    """Create algorithm comparison plots"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_title(f'Algorithm Comparison - {particle_name}', fontsize=14)
    
    # Placeholder bar chart comparing performance metrics
    algorithms = ['PandoraPFA', 'ArborPFA', 'Baseline']
    metrics = ['Energy Res.', 'Efficiency', 'Purity']
    
    x = np.arange(len(algorithms))
    width = 0.25
    
    # Placeholder data - will be replaced with actual metrics
    energy_res = [10, 12, 15]  # %
    efficiency = [92, 90, 88]  # %
    purity = [95, 93, 90]  # %
    
    ax.bar(x - width, energy_res, width, label='Energy Resolution [%]', alpha=0.8)
    ax.bar(x, efficiency, width, label='Efficiency [%]', alpha=0.8)
    ax.bar(x + width, purity, width, label='Purity [%]', alpha=0.8)
    
    ax.set_ylabel('Performance [%]')
    ax.set_xticks(x)
    ax.set_xticklabels(algorithms)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(output_dir) / f"algorithm_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="PFA Single Particle Analysis")
    parser.add_argument("--input-pattern", required=True,
                        help="Glob pattern for input ROOT files")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for plots")
    parser.add_argument("--particle", required=True,
                        help="Particle type (e-, pi-, etc.)")
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Analyzing {args.particle} events...")
    print(f"Input pattern: {args.input_pattern}")
    print(f"Output directory: {args.output_dir}")
    
    try:
        # Load events
        events = load_events(args.input_pattern)
        print(f"Loaded {len(events)} events")
        
        # Extract collections
        data = extract_collections(events)
        
        # Create plots
        print("\nGenerating plots...")
        plot_energy_resolution(data, args.particle, output_dir)
        plot_efficiency(data, args.particle, output_dir)
        plot_algorithm_comparison(data, args.particle, output_dir)
        
        # TODO: Add more analysis:
        # - Momentum resolution
        # - Angular resolution
        # - Detailed efficiency vs eta
        # - PID performance
        
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE")
        print("="*60)
        print(f"Results saved to: {output_dir}")
        print("\nNOTE: This is a placeholder implementation.")
        print("Production version requires:")
        print("  1. Confirmed collection names from EICrecon")
        print("  2. Proper truth matching using EDM4eic links")
        print("  3. Complete metric calculations")
        
    except Exception as e:
        print(f"\nERROR: Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
