#!/bin/bash

source /etc/profile

ddsim --runType batch -N 300 \
	--inputFiles ./data/forward_ions.hepmc \
	--compactFile ./trackers/roman_pot.xml \
	--outputFile ./sim_output/roman_pot_out.root
