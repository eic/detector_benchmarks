#!/bin/bash

ddsim --runType batch --numberOfEvents 10000 \
      --filter.tracker edep0 \
      --compactFile benchmarks/calorimeters/topside.xml \
      --inputFiles  ./data/emcal_pions_upto1GeV_10kevents.hepmc \
      --outputFile  ./sim_output/sim_crystal_pion_input.edm4hep.root
