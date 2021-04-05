#!/bin/bash

if [[ ! -n  "${JUGGLER_DETECTOR}" ]] ; then 
  export JUGGLER_DETECTOR="topside"
fi

if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=1000
fi

if [[ ! -n  "${JUGGLER_INSTALL_PREFIX}" ]] ; then 
  export JUGGLER_INSTALL_PREFIX="/usr/local"
fi

if [[ ! -n  "${E_start}" ]] ; then
  export E_start=5.0
fi

if [[ ! -n  "${E_end}" ]] ; then
  export E_end=5.0
fi

export JUGGLER_FILE_NAME_TAG="emcal_barrel_uniform_pions"
export JUGGLER_GEN_FILE="${JUGGLER_FILE_NAME_TAG}.hepmc"

export JUGGLER_SIM_FILE="sim_${JUGGLER_FILE_NAME_TAG}.root"
export JUGGLER_REC_FILE="rec_${JUGGLER_FILE_NAME_TAG}.root"

echo "JUGGLER_N_EVENTS = ${JUGGLER_N_EVENTS}"
echo "JUGGLER_DETECTOR = ${JUGGLER_DETECTOR}"

# Generate the input events
root -b -q "calorimeters/scripts/emcal_barrel_pions.cxx(${JUGGLER_N_EVENTS}, ${E_start}, ${E_end}, \"${JUGGLER_FILE_NAME_TAG}.hepmc\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: generating input events"
  exit 1
fi
# Plot the input events
root -b -q "calorimeters/scripts/emcal_barrel_pions_reader.cxx(${E_start}, ${E_end}, \"${JUGGLER_FILE_NAME_TAG}.hepmc\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: plotting input events"
  exit 1
fi

# Run geant4 simulations
npsim --runType batch \
      -v WARNING \
      --part.minimalKineticEnergy 0.5*GeV  \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --compactFile topside/${JUGGLER_DETECTOR}.xml \
      --inputFiles ${JUGGLER_FILE_NAME_TAG}.hepmc \
      --outputFile sim_output/${JUGGLER_SIM_FILE}

if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running npdet"
  exit 1
fi

# Run Juggler
xenv -x ${JUGGLER_INSTALL_PREFIX}/Juggler.xenv \
        gaudirun.py calorimeters/options/emcal_barrel_reco.py
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running juggler"
  exit 1
fi

# Directory for plots
mkdir -p results

# Move ROOT output file
mv ${JUGGLER_REC_FILE} sim_output/

