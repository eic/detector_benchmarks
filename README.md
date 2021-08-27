EIC Detector Benchmarks
=======================

[![pipeline status](https://eicweb.phy.anl.gov/EIC/benchmarks/detector_benchmarks/badges/master/pipeline.svg)](https://eicweb.phy.anl.gov/EIC/benchmarks/detector_benchmarks/-/commits/master)

## Overview

Detector benchmarks are meant to test for regressions in individual detector subsystems.
The analysis is meant to avoid a reconstruction step. 
So this precludes using [juggler](https://eicweb.phy.anl.gov/EIC/juggler) for processing the events.

## Documentation

See [common_bench](https://eicweb.phy.anl.gov/EIC/benchmarks/common_bench/).

## Adding new benchmarks

To get an idea of what to do look at an existing benchmark in the 
[`benchmarks` directory](https://eicweb.phy.anl.gov/EIC/benchmarks/detector_benchmarks/-/tree/master/benchmarks).

## Running Locally

### Local development example

Here we setup to use our local build of the `juggler` library.
Note juggler is not needed for `detector_benchmarks` because it is not used but this is the same setup for 
`reconstruction_benchmarks` and `physics_benchmarks`.

First set some environment variables.
```
# export JUGGLER_INSTALL_PREFIX=$HOME/stow/juggler # not needed for detector_benchmarks
export JUGGLER_DETECTOR=athena   # athena is the default
export BEAMLINE_CONFIG=ip6       # ip6 is the default
```


```
git@eicweb.phy.anl.gov:EIC/benchmarks/detector_benchmarks.git && cd detector_benchmarks
git clone https://eicweb.phy.anl.gov/EIC/benchmarks/common_bench.git setup
source setup/bin/env.sh && ./setup/bin/install_common.sh
source .local/bin/env.sh && build_detector.sh
mkdir_local_data_link sim_output
mkdir -p results
mkdir -p config
```




### Pass/Fail tests

- Create a script that returns exit status 0 for success.
- Any non-zero value will be considered failure.
- Script  

