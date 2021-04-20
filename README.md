EIC Detector Benchmarks
=======================

[![pipeline status](https://eicweb.phy.anl.gov/jihee.kim/benchmarks/badges/master/pipeline.svg)](https://eicweb.phy.anl.gov/jihee.kim/benchmarks/-/commits/master)

## Overview

Detector benchmarks are meant to test for regressions in individual detector subsystems.
The analysis is meant to avoid a reconstruction step. 
So this precludes using [juggler](https://eicweb.phy.anl.gov/EIC/juggler) for processing the events.

## Adding new benchmarks

### Pass/Fail tests

- Create a script that returns exit status 0 for success.
- Any non-zero value will be considered failure.
- Script  

