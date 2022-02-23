#!/bin/bash

ddsim --runType batch --numberOfEvents 100 \
      --compactFile benchmarks/crystal_calorimeter/Crystal_example.xml \
      --inputFiles  ./data/emcal_electrons.hepmc \
      --outputFile  ./sim_output/output_emcal_electrons.edm4hep.root
