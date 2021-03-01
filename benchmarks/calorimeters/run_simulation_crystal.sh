#!/bin/bash

ddsim --runType batch --numberOfEvents 100 \
      --compactFile ./calorimeters/Crystal_example.xml \
      --inputFiles  ./data/emcal_electrons.hepmc \
      --outputFile  ./sim_output/output_emcal_electrons.root
