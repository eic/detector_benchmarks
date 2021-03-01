#!/bin/bash

npsim --runType batch --numberOfEvents 100 \
      --compactFile benchmarks/calorimeters/ZDC_example.xml \
      --inputFiles  ./data/zdc_photons.hepmc \
      --outputFile  ./sim_output/output_zdc_photons.root
