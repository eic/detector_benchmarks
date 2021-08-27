#!/bin/bash

source .local/bin/env.sh

if [[ ! -n  "${JUGGLER_DETECTOR}" ]] ; then 
  export JUGGLER_DETECTOR="topside"
fi

if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=100
fi

export FILE_NAME_TAG="forward_protons"
export JUGGLER_GEN_FILE="${LOCAL_DATA_PATH}/${FILE_NAME_TAG}.hepmc"
export JUGGLER_SIM_FILE="${LOCAL_DATA_PATH}/sim_${FILE_NAME_TAG}.root"

echo "JUGGLER_N_EVENTS = ${JUGGLER_N_EVENTS}"
echo "JUGGLER_DETECTOR = ${JUGGLER_DETECTOR}"

# Generate the input events
root -b -q "benchmarks/b0_tracker/scripts/gen_${FILE_NAME_TAG}.cxx(${JUGGLER_N_EVENTS}, \"${FILE_NAME_TAG}.hepmc\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: generating input events"
  exit 1
fi

# Run geant4 simulations
npsim --runType batch \
      -v WARNING \
      --part.minimalKineticEnergy 0.5*GeV  \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --compactFile ${DETECTOR_PATH}/${JUGGLER_DETECTOR}.xml \
      --inputFiles ${FILE_NAME_TAG}.hepmc \
      --outputFile ${JUGGLER_SIM_FILE}

if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running npdet"
  exit 1
fi

# Directory for plots
mkdir -p results

rootls -t ${JUGGLER_SIM_FILE}
# Plot the input events
root -b -q "benchmarks/b0_tracker/scripts/b0_tracker_hits.cxx(\"${JUGGLER_SIM_FILE}\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: events"
  exit 1
fi

