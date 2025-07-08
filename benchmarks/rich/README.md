# RICH Benchmarks

## Setup
- Build `reconstruction_benchmarks` with `cmake`
  - Builds code in `src/` and `include/` to a library
    - Benchmarks are designed as independent algorithms, similar to reconstruction algorithms
    - Algorithm inputs are PODIO collections produced from reconstruction, and the output
      is typically a set of `ROOT` plots
  - Builds executable `benchmark_rich_reconstruction`, which runs any or all of the benchmark
    algorithms and streams any plots to an output `ROOT` file

## Running
- Run the benchmark executable using the wrapper script `run_benchmark.rb`
  - Run with no arguments for usage guide
  - This script can run simulation, reconstruction and/or benchmarks
  - Run `benchmark_rich_reconstruction` to just run the benchmark algorithms
