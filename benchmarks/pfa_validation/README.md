# Particle Flow Analysis (PFA) Validation Benchmark

## Overview

This benchmark validates the dual-track particle flow analysis (PFA) system implemented in EICrecon, comparing:
- **PandoraPFA** - Standard Pandora reconstruction algorithm
- **ArborPFA** - Alternative Arbor-based reconstruction algorithm
- **Baseline** - Existing particle_flow algorithms (for reference)

## Purpose

Comprehensive physics validation of the PFA algorithms through:
1. Single particle validation (energy/momentum resolution, efficiency)
2. Multi-particle validation (multiplicity, jet reconstruction, particle separation)
3. Cross-algorithm performance comparison

## Quick Start

From the top-level detector_benchmarks directory:

```bash
# Run full local validation
snakemake -c2 results/pfa_validation/local

# Or by rule name
snakemake --cores 2 pfa_validation_local

# Process specific campaign
snakemake -c2 results/pfa_validation/24.04.0
```

## Event Types

### Single Particle Events
- **Particles:** electrons (e±), pions (π±), photons (γ)
- **Energies:** 1, 5, 10, 20, 50 GeV
- **Angular regions:**
  - Forward: 2-28° (forward endcap)
  - Barrel: 45-135° (barrel region)
  - Backward: 130-177° (backward endcap)
- **Statistics:** 10,000 events per configuration

### Multi-Particle Events
- **Di-jet events:** Quark-antiquark pairs at various energies
- **Multi-hadron events:** Mixed π±, K±, protons
- **EM showers:** e±, γ combinations
- **Statistics:** 1,000 events per configuration

## Validation Metrics

### Single Particle Validation
1. **Energy Resolution:** (E_reco - E_true) / E_true vs energy
2. **Momentum Resolution:** |p_reco - p_true| / p_true vs momentum
3. **Angular Resolution:** Δθ, Δφ distributions
4. **Efficiency:** Reconstruction efficiency vs energy and η
5. **Particle ID:** PID assignment accuracy

### Multi-Particle Validation
1. **Multiplicity:** Reconstructed vs true particle count
2. **Jet Energy Resolution:** For di-jet events
3. **Missing Energy:** Transverse missing energy resolution
4. **Particle Separation:** Separation efficiency vs ΔR
5. **Event Metrics:** Total visible energy, energy balance

## Output Structure

```
results/pfa_validation/local/
├── single_particle/
│   ├── electrons/
│   │   ├── energy_resolution_vs_energy.png
│   │   ├── efficiency_vs_energy.png
│   │   ├── efficiency_vs_eta.png
│   │   └── algorithm_comparison.png
│   ├── pions/
│   │   └── (similar plots)
│   └── photons/
│       └── (similar plots)
├── multiparticle/
│   ├── dijets/
│   │   ├── jet_energy_resolution.png
│   │   ├── multiplicity_comparison.png
│   │   └── algorithm_performance.png
│   └── hadrons/
│       └── (similar plots)
└── comparison/
    ├── pandora_vs_arbor_energy.png
    ├── pandora_vs_arbor_efficiency.png
    ├── all_algorithms_energy.png
    └── performance_matrix.png
```

## Pass/Fail Criteria

The benchmark passes if:
- ✅ Energy resolution < 15% for E > 5 GeV (single particles)
- ✅ Reconstruction efficiency > 90% for particles in acceptance
- ✅ PandoraPFA and ArborPFA produce valid outputs (no crashes)
- ✅ Multiplicity reconstruction within 20% of truth (multi-particle events)
- ✅ No major regression vs baseline particle_flow algorithms

## Dependencies

- **EICrecon branch:** `wdconinc-pandora-arbor-pfa-with-xml-override`
- **Python packages:** uproot, awkward, matplotlib, numpy, scipy
- **ROOT:** For event generation
- **Simulation:** npsim/ddsim
- **Environment:** eic-shell

## Implementation Details

### Workflow Pipeline
1. **Event Generation:** ROOT macros generate HepMC files with specified kinematics
2. **Simulation:** npsim processes HepMC → EDM4hep (detector simulation)
3. **Reconstruction:** eicrecon processes EDM4hep → EDM4eic (with PFA algorithms enabled)
4. **Analysis:** Python scripts extract collections, compute metrics, generate plots

### Analysis Scripts
- `analysis/gen_particles.cxx` - Single particle event generation
- `analysis/gen_multiparticle.cxx` - Multi-particle event generation
- `analysis/pfa_single_particle.py` - Single particle analysis
- `analysis/pfa_multiparticle.py` - Multi-particle analysis
- `analysis/pfa_comparison.py` - Cross-algorithm comparison

## Parameter Tuning

The PFA algorithms have runtime-configurable parameters. To test different configurations:

1. Modify eicrecon command in Snakefile to include parameter overrides:
   ```bash
   eicrecon input.root -Ppodio:output_file=output.root \
     -PPandoraPFA:SomeParameter=value \
     -PArborPFA:SomeOtherParameter=value
   ```

2. Document parameter changes in validation plots
3. Compare performance across parameter sets

## Troubleshooting

### Common Issues

**Issue:** Reconstruction crashes with "collection not found"
- **Solution:** Verify correct EICrecon branch is checked out
- Check output collection names match configuration

**Issue:** Empty validation plots
- **Solution:** Check reconstruction logs for algorithm failures
- Verify truth matching is working correctly

**Issue:** Low reconstruction efficiency
- **Solution:** Review detector acceptance cuts
- Check if particles are within fiducial volume

**Issue:** Large energy resolution discrepancies
- **Solution:** Verify calibration constants are loaded
- Check if clustering algorithms are properly initialized

## Contributing

When modifying this benchmark:
1. Test full pipeline end-to-end in clean environment
2. Update pass/fail criteria if needed
3. Document any parameter changes
4. Regenerate reference plots
5. Update this README

## References

- EICrecon PFA implementation: https://github.com/eic/EICrecon
- Pandora documentation: https://pandorapfa.github.io/
- Detector benchmarks tutorial: https://eic.github.io/tutorial-developing-benchmarks/

## Authors

- Benchmark implementation: GitHub Copilot (assisted development)
- PFA algorithms: EIC reconstruction team
- Validation framework: detector_benchmarks contributors

## Contact

For questions or issues:
- Open issue on eic/detector_benchmarks repository
- Contact EIC reconstruction team via GitHub
