#!/bin/bash

source /usr/local/bin/thisdd4hep.sh

ddsim --runType batch --numberOfEvents 100 \
      --compactFile ./ZDC_example.xml \
      --inputFiles  ./data/zdc_photons.hepmc \
      --outputFile  ./sim_output/output_zdc_photons.root
