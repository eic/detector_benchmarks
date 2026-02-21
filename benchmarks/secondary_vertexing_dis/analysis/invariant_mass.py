#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description='Analyze invariant mass from secondary vertex data')
    parser.add_argument('--input', '-i', nargs='+', default=['podio_output.root'], 
                       help='Input ROOT file(s) (default: podio_output.root)')
    parser.add_argument('--output', '-o', default='invariant_mass.png',
                       help='Output plot filename (default: invariant_mass.png)')
    parser.add_argument('--vertices-collection', default='SecondaryVerticesHelix',
                       help='Name of the vertices collection (default: SecondaryVerticesHelix)')
    args = parser.parse_args()

    print(f"Loading data from {len(args.input)} file(s): {args.input}")
    
    # Combine data from all input files
    all_sec_x_raw = []
    all_sec_y_raw = []
    all_sec_z_raw = []
    all_sec_chi2_raw = []
    all_sec_ndf_raw = []
    all_px_raw = []
    all_py_raw = []
    all_pz_raw = []
    all_energy_raw = []
    all_assoc_begin = []
    all_assoc_end = []
    all_assoc_indices = []
    
    for input_file in args.input:
        print(f"  Loading {input_file}...")
        file = uproot.open(input_file)
        tree = file['events']
        
        # Load data from this file
        vertices_collection = args.vertices_collection
        all_sec_x_raw.extend(tree[f'{vertices_collection}/{vertices_collection}.position.x'].array())
        all_sec_y_raw.extend(tree[f'{vertices_collection}/{vertices_collection}.position.y'].array())
        all_sec_z_raw.extend(tree[f'{vertices_collection}/{vertices_collection}.position.z'].array())
        all_sec_chi2_raw.extend(tree[f'{vertices_collection}/{vertices_collection}.chi2'].array())
        all_sec_ndf_raw.extend(tree[f'{vertices_collection}/{vertices_collection}.ndf'].array())
        all_px_raw.extend(tree['ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.x'].array())
        all_py_raw.extend(tree['ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.y'].array())
        all_pz_raw.extend(tree['ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.z'].array())
        all_energy_raw.extend(tree['ReconstructedChargedParticles/ReconstructedChargedParticles.energy'].array())
        all_assoc_begin.extend(tree[f'{vertices_collection}/{vertices_collection}.associatedParticles_begin'].array())
        all_assoc_end.extend(tree[f'{vertices_collection}/{vertices_collection}.associatedParticles_end'].array())
        all_assoc_indices.extend(tree[f'_{vertices_collection}_associatedParticles/_{vertices_collection}_associatedParticles.index'].array())
        
        file.close()
    
    print(f"Loaded data from {len(all_sec_x_raw)} events")

    # Calculate total statistics
    total_events = len(all_sec_x_raw)
    total_sec_vertices = sum(len(all_sec_x_raw[i]) for i in range(len(all_sec_x_raw)))
    total_particles = sum(len(all_px_raw[i]) for i in range(len(all_px_raw)))
    
    print(f"\n=== INITIAL STATISTICS ===")
    print(f"Total events: {total_events}")
    print(f"Total secondary vertices: {total_sec_vertices}")
    print(f"Total reconstructed particles: {total_particles}")
    print(f"Total association entries: {sum(len(all_assoc_indices[i]) for i in range(len(all_assoc_indices)))}")

    # Process events with correct associations
    invariant_masses = []
    decay_lengths_2d = []
    decay_lengths_3d = []
    vertex_quality = []
    n_tracks_per_vertex = []

    # Cut flow counters
    events_processed = 0
    vertices_processed = 0
    two_track_vertices = 0
    valid_associations = 0
    valid_invariant_masses = 0
    passed_quality_cuts = 0
    kaon_candidates = 0

    print("\n=== PROCESSING WITH CORRECT ASSOCIATIONS ===")
    
    for event_idx in range(len(all_sec_x_raw)):
        # Get event data
        event_sec_x = all_sec_x_raw[event_idx]
        event_sec_y = all_sec_y_raw[event_idx] 
        event_sec_z = all_sec_z_raw[event_idx]
        event_chi2 = all_sec_chi2_raw[event_idx]
        event_ndf = all_sec_ndf_raw[event_idx]
        
        event_px = all_px_raw[event_idx]
        event_py = all_py_raw[event_idx]
        event_pz = all_pz_raw[event_idx]
        event_energy = all_energy_raw[event_idx]
        
        event_assoc_begin = all_assoc_begin[event_idx]
        event_assoc_end = all_assoc_end[event_idx]
        event_assoc_indices = all_assoc_indices[event_idx]
        
        # Skip if no particles or vertices
        if len(event_px) == 0 or len(event_sec_x) == 0:
            continue
            
        events_processed += 1
            
        # Process each secondary vertex in this event
        for vtx_idx in range(len(event_sec_x)):
            vertices_processed += 1
            
            vtx_begin = event_assoc_begin[vtx_idx]
            vtx_end = event_assoc_end[vtx_idx]
            
            # Get the track indices for this vertex using correct slicing
            if vtx_begin < len(event_assoc_indices) and vtx_end <= len(event_assoc_indices):
                vertex_track_indices = event_assoc_indices[vtx_begin:vtx_end]
                n_tracks = len(vertex_track_indices)
                n_tracks_per_vertex.append(n_tracks)
                
                print(f"Event {event_idx}, Vertex {vtx_idx}: {n_tracks} tracks, indices: {vertex_track_indices}")
                
                # Calculate decay length for this vertex
                decay_2d = np.sqrt(event_sec_x[vtx_idx]**2 + event_sec_y[vtx_idx]**2)
                decay_3d = np.sqrt(event_sec_x[vtx_idx]**2 + event_sec_y[vtx_idx]**2 + event_sec_z[vtx_idx]**2)
                chi2_ndf = event_chi2[vtx_idx] / (event_ndf[vtx_idx] + 1e-10)
                
                # Focus on 2-track vertices (typical for kaon decays)
                if n_tracks == 2:
                    two_track_vertices += 1
                    
                    idx1, idx2 = vertex_track_indices[0], vertex_track_indices[1]
                    
                    # Check if indices are valid for this event
                    if idx1 < len(event_px) and idx2 < len(event_px):
                        valid_associations += 1
                        
                        # Calculate invariant mass
                        p1 = [event_px[idx1], event_py[idx1], event_pz[idx1], event_energy[idx1]]
                        p2 = [event_px[idx2], event_py[idx2], event_pz[idx2], event_energy[idx2]]
                        
                        total_E = p1[3] + p2[3]
                        total_px = p1[0] + p2[0]
                        total_py = p1[1] + p2[1]
                        total_pz = p1[2] + p2[2]
                        
                        inv_mass_sq = total_E**2 - (total_px**2 + total_py**2 + total_pz**2)
                        
                        if inv_mass_sq > 0:
                            valid_invariant_masses += 1
                            inv_mass = np.sqrt(inv_mass_sq)
                            
                            # Apply kaon selection cuts
                            passes_vertex_quality = chi2_ndf < 3.0  # Good vertex fit
                            passes_decay_length = 0.05 < decay_2d < 2.0  # Reasonable decay length for kaons
                            
                            # Only store candidates that pass quality cuts
                            if passes_vertex_quality and passes_decay_length:
                                passed_quality_cuts += 1
                                invariant_masses.append(inv_mass)
                                decay_lengths_2d.append(decay_2d)
                                decay_lengths_3d.append(decay_3d)
                                vertex_quality.append(chi2_ndf)
                                
                                # Check if it's a kaon candidate
                                if 0.4 <= inv_mass <= 0.6:
                                    kaon_candidates += 1
                                    print(f"  -> KAON CANDIDATE: mass={inv_mass:.4f} GeV, decay_length={decay_2d:.3f} mm, χ²/ndf={chi2_ndf:.2f}")
                            else:
                                print(f"  -> Failed cuts: χ²/ndf={chi2_ndf:.2f} (cut<3.0), decay_length={decay_2d:.3f}mm (cut: 0.05-2.0mm)")
                        else:
                            print(f"  -> Invalid mass² = {inv_mass_sq}")
                    else:
                        print(f"  -> Invalid track indices: {idx1}, {idx2} (max: {len(event_px)-1})")
                elif n_tracks > 0:
                    print(f"  -> {n_tracks}-track vertex (not 2-track)")

    # Convert to numpy arrays
    invariant_masses = np.array(invariant_masses)
    decay_lengths_2d = np.array(decay_lengths_2d)
    decay_lengths_3d = np.array(decay_lengths_3d)
    vertex_quality = np.array(vertex_quality)
    n_tracks_per_vertex = np.array(n_tracks_per_vertex)

    print(f"\n=== CUT FLOW RESULTS ===")
    print(f"Events processed:                {events_processed:>8}")
    print(f"Vertices processed:              {vertices_processed:>8}")
    print(f"Two-track vertices:              {two_track_vertices:>8} ({100*two_track_vertices/vertices_processed:.1f}% of vertices)")
    print(f"Valid associations:              {valid_associations:>8} ({100*valid_associations/two_track_vertices:.1f}% of 2-track vertices)")
    print(f"Valid invariant masses:          {valid_invariant_masses:>8} ({100*valid_invariant_masses/valid_associations:.1f}% of valid associations)")
    print(f"Passed quality cuts:             {passed_quality_cuts:>8} ({100*passed_quality_cuts/valid_invariant_masses:.1f}% of valid masses)")
    print(f"Kaon candidates (0.4-0.6 GeV):   {kaon_candidates:>8} ({100*kaon_candidates/passed_quality_cuts:.1f}% of quality candidates)")

    if len(invariant_masses) > 0:
        print(f"\n=== ANALYSIS RESULTS (AFTER CUTS) ===")
        print(f"Applied cuts: χ²/ndf < 3.0, 0.05 < decay_length < 2.0 mm")
        print(f"Mass range: {np.min(invariant_masses):.4f} - {np.max(invariant_masses):.4f} GeV")
        print(f"Mean mass: {np.mean(invariant_masses):.4f} ± {np.std(invariant_masses):.4f} GeV")
        print(f"Decay length range: {np.min(decay_lengths_2d):.3f} - {np.max(decay_lengths_2d):.3f} mm")
        print(f"Mean decay length: {np.mean(decay_lengths_2d):.3f} ± {np.std(decay_lengths_2d):.3f} mm")
        print(f"Vertex quality range: {np.min(vertex_quality):.2f} - {np.max(vertex_quality):.2f}")

        # Track multiplicity analysis
        print(f"\n=== TRACK MULTIPLICITY ===")
        unique_multiplicities, counts = np.unique(n_tracks_per_vertex, return_counts=True)
        for mult, count in zip(unique_multiplicities, counts):
            percentage = 100 * count / len(n_tracks_per_vertex)
            print(f"{mult} tracks per vertex: {count:>8} ({percentage:>5.1f}%)")

        # Decay length cuts analysis
        print(f"\n=== DECAY LENGTH CUTS ANALYSIS ===")
        decay_cuts = {
            'Very Short (<0.1mm)': (0, 0.1),
            'Short (0.1-0.5mm)': (0.1, 0.5),
            'Medium (0.5-1.0mm)': (0.5, 1.0),
            'Long (>1.0mm)': (1.0, np.inf)
        }
        
        for cut_name, (min_cut, max_cut) in decay_cuts.items():
            mask = (decay_lengths_2d >= min_cut) & (decay_lengths_2d < max_cut)
            masses = invariant_masses[mask]
            kaon_mask = (masses >= 0.45) & (masses <= 0.55)
            
            print(f"{cut_name:<25}: {len(masses):>6} pairs, {np.sum(kaon_mask):>3} kaon candidates")

        # Create plot
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))

        # Mass vs decay length (2D histogram)
        h = ax.hist2d(invariant_masses, decay_lengths_2d, 
                     bins=[np.linspace(0, 1, 21), 15], 
                     cmap='viridis', alpha=0.8, range=[[0, 1], None])
        plt.colorbar(h[3], ax=ax, label='Count')
        ax.axvline(0.4976, color='red', linestyle='--', linewidth=2, alpha=0.9, label='K⁰_S')
        ax.axvline(0.4937, color='orange', linestyle=':', linewidth=2, alpha=0.9, label='K±')
        ax.axvline(0.7755, color='purple', linestyle='-.', linewidth=2, alpha=0.9, label='ρ⁰')
        ax.set_xlabel('Invariant Mass (GeV)')
        ax.set_ylabel('2D Decay Length (mm)')
        ax.set_title('Mass vs Decay Length (2D Histogram)\n(0-1 GeV Mass Range)')
        ax.set_xlim(0, 1)
        ax.legend()

        plt.tight_layout()
        plt.savefig(args.output, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {args.output}")

    else:
        print("No valid invariant masses calculated!")

if __name__ == "__main__":
    main()
