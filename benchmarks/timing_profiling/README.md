# Detector Timing Profiling Benchmark

## Overview

This benchmark measures the computational performance of detector simulation by profiling the time spent in different detector regions during Geant4 stepping. It uses a custom `PerformanceProfileSteppingAction` plugin for `npsim` to collect timing data at each simulation step.

## Purpose

- **Identify computational hotspots** in detector geometry
- **Quantify time spent** in different detector subsystems
- **Optimize detector simulation** by focusing on the most expensive regions
- **Track performance** over time as detector designs evolve

## How It Works

1. **Simulation**: Runs `npsim` with the `PerformanceProfileSteppingAction` plugin, which:
   - Records the time spent at each Geant4 step
   - Creates 2D histograms (ZR and XY) of step locations weighted by time
   - Outputs timing data to ROOT histograms

2. **Analysis**: Processes the timing histograms to:
   - Identify the most computationally expensive detector regions
   - Generate annotated plots showing hotspots
   - Create CSV summary of timing by region

## Prerequisites

The `PerformanceProfileSteppingAction` plugin must be built in the `npsim` directory:

```bash
cd npsim
mkdir -p build && cd build
cmake ..
make
```

## Running the Benchmark

### Standalone

```bash
cd detector_benchmarks
bash benchmarks/timing_profiling/run_profile.sh
```

Environment variables (optional):
- `JUGGLER_N_EVENTS`: Number of events to simulate (default: 10)
- `DETECTOR_CONFIG`: Detector configuration to use (default: epic_craterlake)
- `EBEAM`: Electron beam energy in GeV (default: 18)
- `PBEAM`: Proton beam energy in GeV (default: 275)

### In CI/CD Pipeline

The benchmark is defined in `config.yml` with three stages:
1. `sim:timing_profiling` - Runs the profiling simulation
2. `bench:timing_profiling` - Generates plots and analysis
3. `collect_results:timing_profiling` - Collects results

## Output

### Files Generated

- **Histograms**: `sim_output/timing_profile_<timestamp>.histos.root`
  - `m_zr`: 2D histogram (Z vs R) of time per region
  - `m_xy`: 2D histogram (X vs Y) of time per region

- **Plots**: In `results/timing_profiling/`
  - `timing_profile_zr_annotated.{png,pdf}`: Annotated ZR profile showing hotspots
  - `timing_profile_xy.{png,pdf}`: XY profile visualization
  - `timing_profile_regions.csv`: CSV summary of timing by detector region

### Interpretation

The analysis identifies and ranks detector regions by computational cost:
- **Time per event (ms/evt)**: How much CPU time each region consumes per event
- **Percentage of total**: Fraction of total simulation time spent in each region
- Annotated plots highlight the top 5 most expensive regions

## Detector Regions Analyzed

The script analyzes the following detector regions:
- Forward EM+Hadron Calorimeter
- Dual RICH (dRICH)
- DIRC Bar (straight and triangle sections)
- Electron Endcap Electromagnetic Calorimeter (EEEMCal)
- Barrel Imaging and SciFi Calorimeter

## Example Usage

```bash
# Run profiling with 100 events
export JUGGLER_N_EVENTS=100
bash benchmarks/timing_profiling/run_profile.sh

# Analyze the results
python3 benchmarks/timing_profiling/analysis/plot_timing_profile.py \
    sim_output/timing_profile_<timestamp>.histos.root \
    results/timing_profiling/my_analysis \
    100
```

## Notes

- The profiling adds minimal overhead to simulation time
- Larger event samples provide more statistically significant results
- Results can vary based on physics process and beam energies
- This benchmark uses DIS (Deep Inelastic Scattering) events by default
