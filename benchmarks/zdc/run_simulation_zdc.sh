#!/bin/bash

ddsim --runType batch --numberOfEvents 10 \
      --compactFile benchmarks/zdc/ZDC_example.xml \
      --inputFiles  ./data/zdc_photons.hepmc \
      --outputFile  ./sim_output/output_zdc_photons.edm4hep.root
