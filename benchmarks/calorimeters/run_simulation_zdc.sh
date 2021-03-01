#!/bin/bash

ddsim --runType batch --numberOfEvents 100 \
      --compactFile ./calorimeters/ZDC_example.xml \
      --inputFiles  ./data/zdc_photons.hepmc \
      --outputFile  ./sim_output/output_zdc_photons.root
