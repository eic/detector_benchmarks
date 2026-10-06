#!/usr/bin/env python3
"""
PFA Algorithm Comparison Analysis

Creates comprehensive comparison plots across PandoraPFA, ArborPFA, and baseline
particle_flow algorithms, consolidating results from single and multi-particle analyses.
"""

import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json


def load_summary_data(single_dir, multi_dir):
    """
    Load summary statistics from single and multi-particle analyses
    
    TODO: Implement reading of summary JSON files from previous analysis steps
    For now, using placeholder data
    """
    # Placeholder data structure
    data = {
        "single_particle": {
            "electrons": {
                "pandora": {"energy_res": 0.10, "efficiency": 0.92, "purity": 0.95},
                "arbor": {"energy_res": 0.12, "efficiency": 0.90, "purity": 0.93},
                "baseline": {"energy_res": 0.15, "efficiency": 0.88, "purity": 0.90},
            },
            "pions": {
                "pandora": {"energy_res": 0.12, "efficiency": 0.89, "purity": 0.92},
                "arbor": {"energy_res": 0.14, "efficiency": 0.87, "purity": 0.90},
                "baseline": {"energy_res": 0.18, "efficiency": 0.85, "purity": 0.88},
            },
            "photons": {
                "pandora": {"energy_res": 0.08, "efficiency": 0.94, "purity": 0.96},
                "arbor": {"energy_res": 0.10, "efficiency": 0.92, "purity": 0.94},
                "baseline": {"energy_res": 0.13, "efficiency": 0.90, "purity": 0.91},
            },
        },
        "multi_particle": {
            "pandora": {"mult_res": 0.15, "energy_balance": 0.08, "met_res": 0.20},
            "arbor": {"mult_res": 0.18, "energy_balance": 0.10, "met_res": 0.22},
            "baseline": {"mult_res": 0.22, "energy_balance": 0.12, "met_res": 0.25},
        }
    }
    
    return data


def plot_pandora_vs_arbor(data, output_dir):
    """Direct comparison between PandoraPFA and ArborPFA"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('PandoraPFA vs ArborPFA Comparison', fontsize=14)
    
    # Energy resolution comparison
    ax = axes[0, 0]
    particles = ['Electrons', 'Pions', 'Photons']
    pandora_res = [10, 12, 8]  # placeholder %
    arbor_res = [12, 14, 10]
    
    x = np.arange(len(particles))
    width = 0.35
    ax.bar(x - width/2, pandora_res, width, label='PandoraPFA', alpha=0.8, color='blue')
    ax.bar(x + width/2, arbor_res, width, label='ArborPFA', alpha=0.8, color='red')
    ax.set_ylabel('Energy Resolution [%]')
    ax.set_title('Energy Resolution by Particle Type')
    ax.set_xticks(x)
    ax.set_xticklabels(particles)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    
    # Efficiency comparison
    ax = axes[0, 1]
    pandora_eff = [92, 89, 94]
    arbor_eff = [90, 87, 92]
    
    ax.bar(x - width/2, pandora_eff, width, label='PandoraPFA', alpha=0.8, color='blue')
    ax.bar(x + width/2, arbor_eff, width, label='ArborPFA', alpha=0.8, color='red')
    ax.set_ylabel('Efficiency [%]')
    ax.set_title('Reconstruction Efficiency')
    ax.set_xticks(x)
    ax.set_xticklabels(particles)
    ax.axhline(90, color='gray', linestyle='--', alpha=0.5, label='90% threshold')
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    
    # Multi-particle metrics
    ax = axes[1, 0]
    metrics = ['Multiplicity\nRes.', 'Energy\nBalance', 'Missing ET\nRes.']
    pandora_multi = [15, 8, 20]
    arbor_multi = [18, 10, 22]
    
    x = np.arange(len(metrics))
    ax.bar(x - width/2, pandora_multi, width, label='PandoraPFA', alpha=0.8, color='blue')
    ax.bar(x + width/2, arbor_multi, width, label='ArborPFA', alpha=0.8, color='red')
    ax.set_ylabel('Resolution / Error [%]')
    ax.set_title('Multi-Particle Performance')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    
    # Overall score
    ax = axes[1, 1]
    categories = ['Energy\nResolution', 'Efficiency', 'Purity', 'Multi-Particle']
    pandora_scores = [85, 92, 94, 82]
    arbor_scores = [80, 90, 92, 78]
    
    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
    pandora_scores_plot = pandora_scores + [pandora_scores[0]]
    arbor_scores_plot = arbor_scores + [arbor_scores[0]]
    angles_plot = angles + [angles[0]]
    
    ax = plt.subplot(2, 2, 4, projection='polar')
    ax.plot(angles_plot, pandora_scores_plot, 'o-', linewidth=2, label='PandoraPFA', color='blue')
    ax.fill(angles_plot, pandora_scores_plot, alpha=0.25, color='blue')
    ax.plot(angles_plot, arbor_scores_plot, 'o-', linewidth=2, label='ArborPFA', color='red')
    ax.fill(angles_plot, arbor_scores_plot, alpha=0.25, color='red')
    ax.set_xticks(angles)
    ax.set_xticklabels(categories)
    ax.set_ylim(0, 100)
    ax.set_title('Overall Performance Radar', y=1.08)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    ax.grid(True)
    
    plt.tight_layout()
    output_path = Path(output_dir) / "pandora_vs_arbor_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_all_algorithms(data, output_dir):
    """Comparison of all three algorithms"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('All Algorithms Comparison', fontsize=14)
    
    algorithms = ['PandoraPFA', 'ArborPFA', 'Baseline']
    colors = ['blue', 'red', 'green']
    
    # Energy resolution by particle type
    ax = axes[0, 0]
    particles = ['e-', 'π-', 'γ']
    x = np.arange(len(particles))
    width = 0.25
    
    for i, alg in enumerate(algorithms):
        # Placeholder data
        values = [10 + i*2, 12 + i*2, 8 + i*2]
        ax.bar(x + i*width, values, width, label=alg, alpha=0.8, color=colors[i])
    
    ax.set_ylabel('Energy Resolution [%]')
    ax.set_title('Energy Resolution Comparison')
    ax.set_xticks(x + width)
    ax.set_xticklabels(particles)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    
    # Efficiency comparison
    ax = axes[0, 1]
    for i, alg in enumerate(algorithms):
        values = [92 - i*2, 89 - i*2, 94 - i*2]
        ax.bar(x + i*width, values, width, label=alg, alpha=0.8, color=colors[i])
    
    ax.set_ylabel('Efficiency [%]')
    ax.set_title('Reconstruction Efficiency')
    ax.set_xticks(x + width)
    ax.set_xticklabels(particles)
    ax.axhline(90, color='gray', linestyle='--', alpha=0.5)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    
    # Processing time (placeholder)
    ax = axes[1, 0]
    proc_time = [100, 120, 80]  # relative to baseline
    ax.bar(algorithms, proc_time, alpha=0.8, color=colors)
    ax.set_ylabel('Processing Time [relative]')
    ax.set_title('Computational Performance')
    ax.axhline(100, color='gray', linestyle='--', alpha=0.5, label='Baseline')
    ax.grid(True, axis='y', alpha=0.3)
    ax.legend()
    
    # Overall ranking
    ax = axes[1, 1]
    metrics = ['Energy\nRes.', 'Efficiency', 'Purity', 'Speed']
    
    # Normalized scores (higher is better)
    pandora = [90, 92, 94, 80]
    arbor = [85, 90, 92, 83]
    baseline = [70, 88, 88, 100]
    
    x = np.arange(len(metrics))
    width = 0.25
    
    ax.bar(x - width, pandora, width, label='PandoraPFA', alpha=0.8, color='blue')
    ax.bar(x, arbor, width, label='ArborPFA', alpha=0.8, color='red')
    ax.bar(x + width, baseline, width, label='Baseline', alpha=0.8, color='green')
    
    ax.set_ylabel('Performance Score')
    ax.set_title('Overall Metric Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    ax.set_ylim([0, 110])
    
    plt.tight_layout()
    output_path = Path(output_dir) / "all_algorithms_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def create_performance_matrix(data, output_dir):
    """Create a performance summary matrix"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_title('PFA Algorithm Performance Matrix', fontsize=14)
    
    # Performance matrix (normalized scores, 0-100)
    algorithms = ['PandoraPFA', 'ArborPFA', 'Baseline']
    metrics = [
        'e- Energy Res.',
        'π- Energy Res.',
        'γ Energy Res.',
        'Efficiency',
        'Purity',
        'Multiplicity',
        'Speed',
    ]
    
    # Placeholder scores matrix
    matrix = np.array([
        [90, 85, 70],  # e- energy res
        [88, 83, 68],  # pi- energy res
        [92, 87, 72],  # gamma energy res
        [92, 90, 88],  # efficiency
        [94, 92, 88],  # purity
        [85, 82, 78],  # multiplicity
        [80, 83, 100], # speed
    ])
    
    im = ax.imshow(matrix, cmap='RdYlGn', aspect='auto', vmin=60, vmax=100)
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(algorithms)))
    ax.set_yticks(np.arange(len(metrics)))
    ax.set_xticklabels(algorithms)
    ax.set_yticklabels(metrics)
    
    # Rotate the tick labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # Add text annotations
    for i in range(len(metrics)):
        for j in range(len(algorithms)):
            text = ax.text(j, i, f'{matrix[i, j]:.0f}',
                          ha="center", va="center", color="black", fontweight='bold')
    
    ax.set_xlabel('Algorithm')
    ax.set_ylabel('Metric')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Performance Score', rotation=270, labelpad=20)
    
    plt.tight_layout()
    output_path = Path(output_dir) / "performance_matrix.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="PFA Algorithm Comparison")
    parser.add_argument("--single-dir", required=True,
                        help="Directory with single particle results")
    parser.add_argument("--multi-dir", required=True,
                        help="Directory with multi-particle results")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for comparison plots")
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Creating comparison plots...")
    print(f"Single particle dir: {args.single_dir}")
    print(f"Multi-particle dir: {args.multi_dir}")
    print(f"Output directory: {args.output_dir}")
    
    try:
        # Load summary data from previous analyses
        data = load_summary_data(args.single_dir, args.multi_dir)
        
        # Create comparison plots
        print("\nGenerating comparison plots...")
        plot_pandora_vs_arbor(data, output_dir)
        plot_all_algorithms(data, output_dir)
        create_performance_matrix(data, output_dir)
        
        print("\n" + "="*60)
        print("COMPARISON ANALYSIS COMPLETE")
        print("="*60)
        print(f"Results saved to: {output_dir}")
        print("\nNOTE: This is a placeholder implementation using mock data.")
        print("Production version will load actual results from previous analysis steps.")
        
    except Exception as e:
        print(f"\nERROR: Comparison analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
