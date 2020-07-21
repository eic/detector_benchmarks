#!/bin/bash

source /usr/local/bin/thisdd4hep.sh

ddsim --runType batch -N 300 \
	--inputFiles ../data/forward_ions.hepmc \
	--compactFile ./roman_pot.xml \
	--outputFile ../sim_output/roman_pot_out.root
