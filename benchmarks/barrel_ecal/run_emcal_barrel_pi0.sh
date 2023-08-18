#!/bin/bash

if [[ ! -n  "${DETECTOR}" ]] ; then 
  export DETECTOR="topside"
fi

if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=1000
fi

if [[ ! -n  "${E_start}" ]] ; then
  export E_start=5.0
fi

if [[ ! -n  "${E_end}" ]] ; then
  export E_end=5.0
fi

export JUGGLER_FILE_NAME_TAG="emcal_barrel_uniform_pi0"
export JUGGLER_GEN_FILE="${JUGGLER_FILE_NAME_TAG}.hepmc"

export JUGGLER_SIM_FILE="sim_${JUGGLER_FILE_NAME_TAG}.edm4hep.root"
export JUGGLER_REC_FILE="rec_${JUGGLER_FILE_NAME_TAG}.root"

echo "JUGGLER_N_EVENTS = ${JUGGLER_N_EVENTS}"
echo "DETECTOR = ${DETECTOR}"

# Generate the input events
root -b -q "benchmarks/barrel_ecal/scripts/emcal_barrel_pi0.cxx(${JUGGLER_N_EVENTS}, ${E_start}, ${E_end}, \"${JUGGLER_FILE_NAME_TAG}.hepmc\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: generating input events"
  exit 1
fi
# Plot the input events
root -b -q "benchmarks/barrel_ecal/scripts/emcal_barrel_pi0_reader.cxx(${E_start}, ${E_end}, \"${JUGGLER_FILE_NAME_TAG}.hepmc\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: plotting input events"
  exit 1
fi

# Run geant4 simulations
ddsim --runType batch \
      -v WARNING \
      --part.minimalKineticEnergy 0.5*GeV  \
      --filter.tracker edep0 \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
      --inputFiles ${JUGGLER_FILE_NAME_TAG}.hepmc \
      --outputFile sim_output/${JUGGLER_SIM_FILE}

if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running npdet"
  exit 1
fi

# Directory for plots
mkdir -p results

# Move ROOT output file
#mv ${JUGGLER_REC_FILE} sim_output/

