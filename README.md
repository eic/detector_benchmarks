ePIC Detector Benchmarks
========================

[![pipeline status](https://eicweb.phy.anl.gov/EIC/benchmarks/detector_benchmarks/badges/master/pipeline.svg)](https://eicweb.phy.anl.gov/EIC/benchmarks/detector_benchmarks/-/commits/master)

## Overview

Detector benchmarks are meant to provide a maintained set of performance plots for individual detector subsystems.

## Documentation

 - See [tutorial](https://eic.github.io/tutorial-developing-benchmarks/)
 - See [common_bench](https://eicweb.phy.anl.gov/EIC/benchmarks/common_bench/).

## Adding new benchmarks

To get an idea of what to do look at an existing benchmark in the 
[`benchmarks` directory](https://github.com/eic/detector_benchmarks/tree/master/benchmarks).
Currently a good reference for Snakemake instrumentation is available in the `tracking\_performances` benchmark.
It relies on single particle simulations that can be either produced on eicweb or downloaded from official campagins.

### File organization

For a minimal benchmark you'll need to add
`benchmarks/<benchmark_name_here>/config.yml` and
`benchmarks/<benchmark_name_here>/Snakemake`, plus the analysis script/macro.
The `Snakefile` has to be included in the root `./Snakefile` of the repository.
That common entry point is needed to ensure that common simulation samples can
be defined to be re-used by several benchmarks at a time.
The `config.yml` will require an include from the `./.gitlab-ci.yml`.

### Pass/Fail tests

 - Create a script that returns exit status 0 for success.
 - Any non-zero value will be considered failure.
