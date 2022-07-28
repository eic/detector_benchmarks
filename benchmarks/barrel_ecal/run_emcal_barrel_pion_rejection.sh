#!/bin/bash

if [[ ! -n  "${DETECTOR}" ]] ; then 
  export DETECTOR="athena"
fi

if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=100
fi

if [[ ! -n  "${E_START}" ]] ; then
  export E_START=1.0
fi

if [[ ! -n  "${E_END}" ]] ; then
  export E_END=18.0
fi

export PARTICLE=$1
if [[ ! -n  "${PARTICLE}" ]] ; then
  export PARTICLE="electron"
fi

export JUGGLER_FILE_NAME_TAG="emcal_barrel_piRej_${PARTICLE}"
export JUGGLER_GEN_FILE="${JUGGLER_FILE_NAME_TAG}.hepmc"

export JUGGLER_SIM_FILE="sim_${JUGGLER_FILE_NAME_TAG}.edm4hep.root"
export JUGGLER_REC_FILE="rec_${JUGGLER_FILE_NAME_TAG}.root"

echo "JUGGLER_N_EVENTS = ${JUGGLER_N_EVENTS}"
echo "DETECTOR = ${DETECTOR}"

# Generate the input events
root -b -q "benchmarks/barrel_ecal/scripts/emcal_barrel_particles_gen.cxx+(${JUGGLER_N_EVENTS}, ${E_START}, ${E_END}, \"${PARTICLE}\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: generating input events"
  exit 1
fi
# Plot the input events
root -b -q "benchmarks/barrel_ecal/scripts/emcal_barrel_particles_reader.cxx+(\"${PARTICLE}\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: plotting input events"
  exit 1
fi

ls -ltRhL

ddsim --runType batch \
      -v WARNING \
      --part.minimalKineticEnergy 0.5*GeV  \
      --filter.tracker edep0 \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
      --inputFiles data/${JUGGLER_FILE_NAME_TAG}.hepmc \
      --outputFile sim_output/${JUGGLER_SIM_FILE}

if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running npdet"
  exit 1
fi

# Directory for plots
mkdir -p results

# Move ROOT output file
#mv ${JUGGLER_REC_FILE} sim_output/

