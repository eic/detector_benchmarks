# Automated Benchmark Plot Generation

Simple system to generate PNG plots from detector benchmarks using Snakemake + Slurm.

## Key Features

- **1 Slurm job per benchmark** - Entire pipeline (generation → simulation → reconstruction → analysis) runs in a single job
- **PNG outputs only** - Auto-converts PDFs to PNG
- **Uses eic-shell** - Runs in container with full EIC software stack
- **No refactoring needed** - Works with existing benchmark code

## Quick Start

### Local Testing

```bash
# In eic-shell
eic-shell
cd detector_benchmarks

# Test one benchmark
snakemake test_benchmark --config benchmark=zdc_neutron

# Run all configured benchmarks (4 cores)
snakemake -j 4 all_plots
```

### Cluster Execution (JLab ifarm)

```bash
ssh ifarm.jlab.org
cd detector_benchmarks
mkdir -p logs/slurm

# Submit to Slurm - each benchmark becomes 1 job
snakemake --profile profiles/slurm all_plots

# Monitor jobs
squeue -u $USER

# Check logs
tail -f logs/slurm/*.out
```

Plots are written directly to `results/` as each job completes.

## Configuration

### Add Benchmarks

Edit `rules/all_plots.smk` and add to `BENCHMARK_PLOTS`:

```python
BENCHMARK_PLOTS = {
    "zdc_neutron": lambda cfg: [
        f"results/zdc_neutron/{cfg}/fwd_neutrons_geant.png",
        f"results/zdc_neutron/{cfg}/fwd_neutrons_recon.png",
    ],
    "your_benchmark": lambda cfg: [
        f"results/your_benchmark/{cfg}/plot.png",
    ],
}
```

### Cluster Settings

Edit `profiles/slurm/config.yaml`:

```yaml
slurm_partition: "production"    # Your partition
slurm_account: null              # Set to "eic" if required
runtime: 240                     # Minutes per job
mem_mb: 4000                     # Memory per job
jobs: 50                         # Max concurrent jobs
```

### Group Benchmark Steps (1 job per benchmark)

To run all steps of a benchmark in a single Slurm job, add `group:` to each rule in the benchmark's Snakefile:

```python
# benchmarks/your_benchmark/Snakefile

rule step1:
    group: "your_benchmark_pipeline"
    ...

rule step2:
    group: "your_benchmark_pipeline"
    ...

rule step3:
    group: "your_benchmark_pipeline"
    ...
```

See `benchmarks/zdc_neutron/Snakefile` for example.

## Collect Plots

```bash
# Gather all PNG files into one directory
python scripts/consolidate_plots.py

# Output in consolidated_plots/
ls consolidated_plots/
```

## File Structure

```
detector_benchmarks/
├── Snakefile                   # Main workflow
├── rules/all_plots.smk        # Plot collection rules
├── profiles/slurm/config.yaml # Slurm settings
├── benchmarks/
│   └── zdc_neutron/           # Example grouped benchmark
└── results/                   # PNG outputs appear here
```

## How It Works

1. **Snakemake** builds dependency graph of all plots
2. **Grouped benchmarks** run entirely in 1 Slurm job:
   - Generate events
   - Run simulation  
   - Run reconstruction
   - Generate analysis plots
   - Convert PDFs → PNG
3. **Multiple benchmarks** run in parallel (up to `jobs` limit)
4. **Output PNG files** are written directly to `results/` by each job

## Troubleshooting

**Commands not found:**
```bash
eic-shell -- which npsim eicrecon root convert
```

**Still getting multiple jobs per benchmark:**
- Check all rules in benchmark have same `group:` name
- Verify `profiles/slurm/config.yaml` has `group-components: {group: 1}`

**Job timeout:**
- Increase `runtime` in `profiles/slurm/config.yaml`

**Need more memory:**
```python
rule expensive_step:
    group: "pipeline"
    resources:
        mem_mb: 16000  # Applies to whole group
```

## Example

```bash
# Test zdc_neutron locally
snakemake test_benchmark --config benchmark=zdc_neutron

# Output written directly by the job:
# results/zdc_neutron/epic_craterlake/fwd_neutrons_geant.png
# results/zdc_neutron/epic_craterlake/fwd_neutrons_recon.png
```

On Slurm: **1 job** runs all 4 steps (hepmc → sim → reco → analysis)
