#!/bin/bash

ddsim --runType batch -N 300 \
	--inputFiles ./data/forward_ions.hepmc \
	--compactFile benchmarks/trackers/roman_pot.xml \
	--outputFile ./sim_output/roman_pot_out.edm4hep.root
