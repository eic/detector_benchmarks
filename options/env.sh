#!/bin/bash

## =============================================================================
## Global configuration variables for the benchmark scripts
## The script defines the following environment variables that are meant to
## be overriden by the Gitlab continuous integration (CI)
##
##  - JUGGLER_DETECTOR:       detector package to be used for the benchmark
##  - JUGGLER_N_EVENTS:       #events processed by simulation/reconstruction
##  - JUGGLER_INSTALL_PREFIX: location where Juggler (digi/recon) is installed
##  - JUGGLER_N_THREADS:      Number of threads/processes to spawn in parallel
##  - JUGGLER_RNG_SEED:       Random seed for the RNG
##
## It also defines the following additional variables for internal usage
##  - LOCAL_PREFIX:           prefix for packages installed during the benchmark
##  - DETECTOR_PREFIX:        prefix for the detector definitions
##  - DETECTOR_PATH:          actual path with the detector definitions
##
## Finally, it makes sure LOCAL_PREFIX and JUGGLER_PREFIX are added to PATH
## and LD_LIBRARY_PATH
## =============================================================================

echo "Setting up the Physics Benchmarks environment"

## =============================================================================
## Default variable definitions, normally these should be set
## by the CI. In case of local development you may want to change these
## in case you would like to modify the detector package or
## number of events to be analyzed during the benchmark

## Detector package to be used during the benchmark process
if [ ! -n  "${JUGGLER_DETECTOR}" ] ; then 
  export JUGGLER_DETECTOR="topside"
fi

if [ ! -n  "${JUGGLER_DETECTOR_VERSION}" ] ; then 
  export JUGGLER_DETECTOR_VERSION="master"
fi


## Number of events that will be processed by the reconstruction
if [ ! -n  "${JUGGLER_N_EVENTS}" ] ; then 
  export JUGGLER_N_EVENTS=100
fi

## Maximum number of threads or processes a single pipeline should use
## (this is not enforced, but the different pipeline scripts should use
##  this to guide the number of parallel processes or threads they 
##  spawn).
if [ ! -n "${JUGGLER_N_THREADS}" ]; then
  export JUGGLER_N_THREADS=10
fi

## Random seed for event generation, should typically not be changed for
## reproductability.
if [ ! -n "${JUGGLER_RNG_SEED}" ]; then
  export JUGGLER_RNG_SEED=1
fi

## Install prefix for juggler, needed to locate the Juggler xenv files.
## Also used by the CI as install prefix for other packages where needed.
## You should not have to touch this. Note that for local usage a different 
## prefix structure is automatically used.
if [ ! -n  "${JUGGLER_INSTALL_PREFIX}" ] ; then 
  export JUGGLER_INSTALL_PREFIX="/usr/local"
fi
## Ensure the juggler prefix is an absolute path
export JUGGLER_INSTALL_PREFIX=`realpath ${JUGGLER_INSTALL_PREFIX}`


## Location of local data for pass data from job to job within pipeline.
## Not saved as artifacts.
if [ ! -n  "${LOCAL_DATA_PATH}" ] ; then 
  export LOCAL_DATA_PATH="/scratch/${CI_PROJECT_NAME}_${CI_PIPELINE_ID}"
fi

## =============================================================================
## Other utility variables that govern how some of the dependent packages
## are built and installed. You should not have to change these.

## local prefix to be used for local storage of packages
## downloaded/installed during the benchmark process
LOCAL_PREFIX=".local"
mkdir -p ${LOCAL_PREFIX}
export LOCAL_PREFIX=`realpath ${LOCAL_PREFIX}`

## detector prefix: prefix for the detector definitions
export DETECTOR_PREFIX="${LOCAL_PREFIX}/detector"
mkdir -p ${DETECTOR_PREFIX}

## detector path: actual detector definition path
export DETECTOR_PATH="${DETECTOR_PREFIX}/${JUGGLER_DETECTOR}"

## build dir for ROOT to put its binaries etc.
export ROOT_BUILD_DIR=$LOCAL_PREFIX/root_build

echo "JUGGLER_DETECTOR:           ${JUGGLER_DETECTOR}"
echo "JUGGLER_DETECTOR_VERSION:   ${JUGGLER_DETECTOR_VERSION}"
echo "JUGGLER_N_EVENTS:           ${JUGGLER_N_EVENTS}"
echo "JUGGLER_N_THREADS:          ${JUGGLER_N_THREADS}"
echo "JUGGLER_RNG_SEED:           ${JUGGLER_RNG_SEED}"
echo "JUGGLER_INSTALL_PREFIX:     ${JUGGLER_INSTALL_PREFIX}"
echo "LOCAL_DATA_PATH:            ${LOCAL_DATA_PATH}"
echo "LOCAL_PREFIX:               ${LOCAL_PREFIX}"
echo "DETECTOR_PREFIX:            ${DETECTOR_PREFIX}"
echo "DETECTOR_PATH:              ${DETECTOR_PATH}"
echo "ROOT_BUILD_DIR:             ${ROOT_BUILD_DIR}"

## =============================================================================
## Setup PATH and LD_LIBRARY_PATH to include our prefixes
echo "Adding JUGGLER_INSTALL_PREFIX and LOCAL_PREFIX to PATH and LD_LIBRARY_PATH"
export PATH=${JUGGLER_INSTALL_PREFIX}/bin:${LOCAL_PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${JUGGLER_INSTALL_PREFIX}/lib:${LOCAL_PREFIX}/lib:${LD_LIBRARY_PATH}

## =============================================================================
## That's all!
echo "Environment setup complete."
