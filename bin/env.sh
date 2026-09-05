#!/bin/bash

## =============================================================================
## Global configuration variables for the benchmark scripts
## The script defines the following environment variables that are meant to
## be overriden by the Gitlab continuous integration (CI)
##
##  - DETECTOR:                detector package to be used for the benchmark
##  - DETECTOR_VERSION:        detector package to be used for the benchmark
##  - BENCHMARK_N_THREADS:     number of threads/processes to spawn in parallel
##  - BENCHMARK_RNG_SEED:      random seed for the RNG
##
## Note: DETECTOR_CONFIG and BENCHMARK_N_EVENTS are managed by Snakemake
##       via snakemake.yml and can be overridden with --config on the CLI.
##
## It also defines the following additional variables for internally usage
##  - LOCAL_PREFIX:           prefix for packages installed during the benchmark
##  - LOCAL_DATA_PATH:        local storage for pipeline jobs
##  - DETECTOR_PATH:          root path to locally installed detector definition xml files
##
## Finally, it makes sure LOCAL_PREFIX is added to PATH and LD_LIBRARY_PATH
## =============================================================================

echo "Setting up the Physics Benchmarks environment"

## =============================================================================
## Default variable definitions, normally these should be set
## by the CI. In case of local development you may want to change these
## in case you would like to modify the detector package or
## number of events to be analyzed during the benchmark

## Detector package to be used during the benchmark process
## If you didn't define DETECTOR already, you have bigger problesm

if [ ! -n  "${DETECTOR}" ] ; then
  echo "ERROR: No DETECTOR defined!" 
  echo "       There is no assumed default detector."
  echo "       Set the environment variable DETECTOR accordingly."
  #export DETECTOR="epic"
fi

# main is the new master
if [ ! -n  "${DETECTOR_VERSION}" ] ; then 
  export DETECTOR_VERSION="main"
fi


## Maximum number of threads or processes a single pipeline should use
## (this is not enforced, but the different pipeline scripts should use
##  this to guide the number of parallel processes or threads they 
##  spawn).
if [ ! -n "${BENCHMARK_N_THREADS}" ]; then
  export BENCHMARK_N_THREADS=10
fi
export ROOT_MAX_THREADS=${BENCHMARK_N_THREADS}

## Random seed for event generation, should typically not be changed for
## reproductability.
if [ ! -n "${BENCHMARK_RNG_SEED}" ]; then
  export BENCHMARK_RNG_SEED=1
fi

## Location of local data for pass data from job to job within pipeline.
## Not saved as artifacts.
## Local /scratch directory is presumed to be writable. 
if [ ! -n  "${LOCAL_DATA_PATH}" ] ; then
  if [ -w /scratch ] ; then
    export LOCAL_DATA_PATH="/scratch/${CI_PROJECT_NAME}_${CI_PIPELINE_ID}"
  else
    echo "/scratch not writable; using $PWD/scratch"
    export LOCAL_DATA_PATH="$PWD/scratch/${CI_PROJECT_NAME}_${CI_PIPELINE_ID}"
  fi
fi
mkdir -p "${LOCAL_DATA_PATH}"
if [ ! -d "${LOCAL_DATA_PATH}" ]; then 
  echo "LOCAL_DATA_PATH (${LOCAL_DATA_PATH}) does not exist!!"
  echo "Creating LOCAL_DATA_PATH=$(pwd)/local_data "
  export LOCAL_DATA_PATH="$(pwd)/local_data"
  mkdir -p "${LOCAL_DATA_PATH}"
fi

## =============================================================================
## Other utility variables that govern how some of the dependent packages
## are built and installed. You should not have to change these.

## local prefix to be used for local storage of packages
## downloaded/installed during the benchmark process
LOCAL_PREFIX=".local"
mkdir -p "${LOCAL_PREFIX}"
export LOCAL_PREFIX=`realpath ${LOCAL_PREFIX}`

## detector path: root path to locally installed detector definition xml files
export DETECTOR_PATH="${LOCAL_PREFIX}/share/${DETECTOR}"

## build dir for ROOT to put its binaries etc.
export ROOT_BUILD_DIR=$LOCAL_PREFIX/root_build

export ROOT_INCLUDE_PATH=${LOCAL_PREFIX}/include:${ROOT_INCLUDE_PATH}

## =============================================================================
## Setup PATH and LD_LIBRARY_PATH to include our prefixes
echo "Adding LOCAL_PREFIX to PATH and LD_LIBRARY_PATH"
export PATH=${LOCAL_PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${LOCAL_PREFIX}/lib:${LD_LIBRARY_PATH}

# Local field maps
mkdir -p ${LOCAL_DATA_PATH}/fieldmaps
ln -sf ${LOCAL_DATA_PATH}/fieldmaps

## =============================================================================
## That's all!
echo "Environment setup complete."
