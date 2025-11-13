#!/bin/bash

#set -o errexit
#set -o pipefail

function print_the_help {
  echo "USAGE: ${0}  [--sim-only]  "
  echo "OPTIONS: "
  echo "  --sim-only    Only run up to the simulation "
  echo "  --analysis    Only run the analysis scripts "
  exit 
}

ANALYSIS_ONLY=
SIM_ONLY=
POSITIONAL=()

while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    -h|--help)
      shift # past argument
      print_the_help
      ;;
    --sim-only)
      shift # past argument
      SIM_ONLY=1
      ;;
    --analysis)
      shift # past argument
      ANALYSIS_ONLY=1
      ;;
    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $1"
      print_the_help
      shift # past argument
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

print_env.sh

## To run the reconstruction, we need the following global variables:
## - DETECTOR:         the detector package we want to use for this benchmark
## - DETECTOR_VERSION: the detector package we want to use for this benchmark
## - DETECTOR_PATH:            full path to the detector definitions
##
## You can ready options/env.sh for more in-depth explanations of the variables
## and how they can be controlled.

if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=100
fi

export JUGGLER_FILE_NAME_TAG="track_hits"
export JUGGLER_GEN_FILE="${LOCAL_DATA_PATH}/${JUGGLER_FILE_NAME_TAG}.hepmc"

export JUGGLER_SIM_FILE="${LOCAL_DATA_PATH}/sim_${JUGGLER_FILE_NAME_TAG}.edm4hep.root"

echo "JUGGLER_N_EVENTS = ${JUGGLER_N_EVENTS}"
echo "DETECTOR = ${DETECTOR}"


if [ -z "${ANALYSIS_ONLY}" ] ; then

  echo "Generating Events"
  ## generate the input events
  root -b -q "benchmarks/tracking_detectors/scripts/gen_track_hits.cxx(${JUGGLER_N_EVENTS}, \"${JUGGLER_FILE_NAME_TAG}.hepmc\")"
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running script"
    exit 1
  fi

  echo "Running geant4 simulation"
  ## run geant4 simulations
  npsim --runType batch \
    --part.minimalKineticEnergy 1000*GeV  \
    --filter.tracker edep0 \
    -v WARNING \
    --numberOfEvents ${JUGGLER_N_EVENTS} \
    --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
    --inputFiles  ${JUGGLER_FILE_NAME_TAG}.hepmc \
    --outputFile  ${JUGGLER_SIM_FILE}
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running script"
    exit 1
  fi

fi


if [ -z "${SIM_ONLY}" ] ; then

  echo "Running analysis scripts"

  mkdir -p results/tracking_detectors
  rootls -t ${JUGGLER_SIM_FILE}
  root -b -q "benchmarks/tracking_detectors/analysis/sim_track_hits.cxx+(\"${JUGGLER_SIM_FILE}\")"
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running root script"
    exit 1
  fi

fi

root_filesize=$(stat --format=%s "${JUGGLER_SIM_FILE}")
if [[ "${JUGGLER_N_EVENTS}" -lt "500" ]] ; then 
  # file must be less than 10 MB to upload
  if [[ "${root_filesize}" -lt "10000000" ]] ; then 
    cp ${JUGGLER_SIM_FILE} results/.
  fi
fi


