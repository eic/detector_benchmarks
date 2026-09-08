# RICH IRT benchmark

Benchmark of the RICH Indirect Ray Tracing (IRT) reconstruction.

Single particles (`pi-`, `kaon-`, `proton`) are shot into the acceptance of one
of the RICH detectors, reconstructed with `eicrecon`, and the resulting IRT
collections are analyzed:

* `{RICH}Tracks` — track momentum and pseudorapidity at the RICH,
* `{RICH}IrtRadiatorInfo` — per-radiator Cherenkov angle, number of
  photoelectrons and number of hits, and the Cherenkov angle versus momentum,
* `{RICH}IrtParticles` vs. `MCParticles` — PID confusion matrix.

## Layout

| Path | Description |
| --- | --- |
| `src/`, `include/` | C++ analysis program (`irt_analysis`), built with CMake |
| `scripts/plot_irt.py` | Produces the PNG plots from the histogram files |
| `Snakefile` | Simulation, reconstruction, analysis and plotting rules |
| `config.yml` | GitLab CI jobs |

## Running

```console
snakemake --cores 1 results/epic_craterlake/irt/DRICH
```

The phase space (polar angle and momentum range), the particle species and the
number of events are configured at the top of the `Snakefile` in
`IRT_SETTINGS`, `IRT_PARTICLES` and `IRT_N_EVENTS`.

The analysis program can also be run standalone:

```console
irt_analysis <DRICH|PFRICH> <rec.edm4eic.root> <hists.root>
```

## Status

This is the first iteration of the benchmark. The names of the IRT collections
produced by `eicrecon` have not been verified against a released version yet,
and the pfRICH chain is not exercised in CI, so the CI jobs are marked
`allow_failure: true`.
