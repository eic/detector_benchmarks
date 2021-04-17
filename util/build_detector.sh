#!/bin/bash

## =============================================================================
## Build and install the JUGGLER_DETECTOR detector package into our local prefix
## =============================================================================

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/..
pushd ${PROJECT_ROOT}

## =============================================================================
## Load the environment variables. To build the detector we need the following
## variables:
##
## - JUGGLER_DETECTOR: the detector package we want to use for this benchmark
## - LOCAL_PREFIX:     location where local packages should be installed
## - DETECTOR_PREFIX:  prefix for the detector definitions 
## - DETECTOR_PATH:    full path for the detector definitions
##                     this is the same as ${DETECTOR_PREFIX}/${JUGGLER_DETECTOR}
##
## You can read options/env.sh for more in-depth explanations of the variables
## and how they can be controlled.
source options/env.sh

## =============================================================================
## Step 1: download/update the detector definitions (if needed)
pushd ${DETECTOR_PREFIX}

## We need an up-to-date copy of the detector
if [ ! -d ${JUGGLER_DETECTOR} ]; then
  echo "Fetching ${JUGGLER_DETECTOR}"
  git clone -b ${JUGGLER_DETECTOR_VERSION} https://eicweb.phy.anl.gov/EIC/detectors/${JUGGLER_DETECTOR}.git
else
  echo "Updating ${JUGGLER_DETECTOR}"
  pushd ${JUGGLER_DETECTOR}
  git pull --ff-only
  popd
fi
## We also need an up-to-date copy of the accelerator. For now this is done
## manually. Down the road we could maybe automize this with cmake
if [ ! -d accelerator ]; then
  echo "Fetching accelerator"
  git clone https://eicweb.phy.anl.gov/EIC/detectors/accelerator.git
else
  echo "Updating accelerator"
  pushd accelerator
  git pull --ff-only
  popd
fi
## Now symlink the accelerator definition into the detector definition
echo "Linking accelerator definition into detector definition"
ln -s -f ${DETECTOR_PREFIX}/accelerator/eic ${DETECTOR_PATH}/eic

## =============================================================================
## Step 2: Compile and install the detector definition
echo "Building and installing the ${JUGGLER_DETECTOR} package"

mkdir -p ${DETECTOR_PREFIX}/build
pushd ${DETECTOR_PREFIX}/build
cmake ${DETECTOR_PATH} -DCMAKE_INSTALL_PREFIX=${LOCAL_PREFIX} && make -j30 install

## =============================================================================
## Step 3: That's all!
echo "Detector build/install complete!"
