#!/bin/bash

ddsim --runType batch -N 100 \
	--inputFiles ../data/forward_ions.hepmc \
	--compactFile ./roman_pot.xml \
	--outputFile ./roman_pot_out.root
