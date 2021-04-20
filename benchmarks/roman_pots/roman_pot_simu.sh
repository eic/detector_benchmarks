#!/bin/bash

if [[ ! -n  "${JUGGLER_DETECTOR}" ]] ; then 
  export JUGGLER_DETECTOR="topside"
fi

if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=1000
fi

ddsim --runType batch -N 300 \

npsim --runType batch \
      -v WARNING \
      --part.minimalKineticEnergy 0.5*GeV  \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --compactFile ${DETECTOR_PATH}/${JUGGLER_DETECTOR}.xml \
	--inputFiles ./data/forward_ions.hepmc \
      --outputFile sim_output/${JUGGLER_SIM_FILE}
