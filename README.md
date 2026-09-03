ePIC Detector Benchmarks
========================

[![pipeline status](https://eicweb.phy.anl.gov/EIC/benchmarks/detector_benchmarks/badges/master/pipeline.svg)](https://eicweb.phy.anl.gov/EIC/benchmarks/detector_benchmarks/-/commits/master)

## Overview

Detector benchmarks are meant to provide a maintained set of performance plots for individual detector subsystems.

## Documentation

 - See [tutorial](https://eic.github.io/tutorial-developing-benchmarks/).

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

### Caching (`cache: True`) and rule `params`

When using Snakemake's `cache: True`, be aware that the cache key (the
"provenance hash") is computed from:

 - the rule's **unformatted** shell command template — i.e. the literal text
   `{wildcards.FOO}`, not the actual substituted value,
 - the job's `params`,
 - the content hash of input files that are *not* produced by another job,
 - the (recursively computed) provenance hash of upstream jobs, for inputs
   that *are* produced by another job.

As a consequence, wildcards that are only used in the `shell:` command or in
`input:`/`output:` paths, and are not otherwise reflected in `params` or in
the content of a non-generated input file, do **not** affect the cache key.
This can cause jobs with different wildcard values (e.g. different
particles/energies) to incorrectly collide onto the same cache entry.

To avoid this, any wildcard that affects the result of a cached rule (or of
an upstream rule it depends on) must be listed explicitly in `params`, e.g.
`PARTICLE = lambda wildcards: wildcards.PARTICLE`.
