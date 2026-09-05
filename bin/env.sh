#!/bin/bash

## =============================================================================
## Environment setup for benchmark scripts
##
## This script defines the following environment variables:
##  - LOCAL_PREFIX:           prefix for packages installed during the benchmark
##  - LOCAL_DATA_PATH:        local storage for pipeline jobs
##  - ROOT_BUILD_DIR:         build dir for ROOT binaries
##  - ROOT_INCLUDE_PATH:      C++ include path for ROOT
##  - PATH and LD_LIBRARY_PATH: updated to include LOCAL_PREFIX directories
## =============================================================================

echo "Setting up the Physics Benchmarks environment"

## =============================================================================
## Default variable definitions

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
