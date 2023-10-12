#!/bin/bash
##########################################################
# Stand-alone cript to generate ACTS material map
# For ACTS v21
# Instructions: https://indico.bnl.gov/event/20842/contributions/81873/attachments/50367/86133/Auto%20Script%20for%20ACTS%20Material%20Map-2.pdf
# Please read through code and comments, then follow the option 1,2,3,4 
# to generate a material map and validation plots.
#
# Shujie Li, Oct 2023
##########################################################
echo "=========FOR EPIC Craterlake========="
# source /opt/detector/setup.sh


## install ACTS
# cmake -B build -S . -DCMAKE_INSTALL_PREFIX=$ATHENA_PREFIX -DCMAKE_CXX_STANDARD=17 -DACTS_BUILD_PLUGIN_DD4HEP=ON -DACTS_BUILD_PLUGIN_JSON=ON -DACTS_BUILD_PLUGIN_TGEO=ON -DACTS_BUILD_EXAMPLES=ON -DACTS_BUILD_EXAMPLES_DD4HEP=ON -DACTS_BUILD_EXAMPLES_GEANT4=ON -DACTS_BUILD_INTEGRATIONTESTS=ON -DACTS_BUILD_UNITTESTS=ON

## in container: 
# /usr/local/bin/this_acts.sh
#
export ACTS_PATH=${PWD}/scripts ## path to example scripts downloaded from acts/Examples/Scripts/MaterialMapping
export ACTS_BUILD_DIR="/usr/local/"

## the xml file should include all detectors and materials 
## within the volume of interest. Objects outside of this
## volume (e.g. solenoid) can be removed for faster simulation.

# export DETECTOR=epic
export DETECTOR_XML_PATH=${PWD}  #$DETECTOR_PATH
export DETECTOR_NAME="epic_craterlake_no_beampipe"

nev=1000  # number of events for geantino scan
nparticles=5000 # number of particles per event

kopt=$1
while [[ $# -gt 1 ]]
do
  key="$2"

  case $key in
    --nev)
      nev=$3
      shift # past value
      shift
      ;;
    --nparticles)
      nparticles=$3
      shift # past value
      shift
      ;;
    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $2"
      print_the_help
      shift # past argument
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters



## --------------------------
## GEANTINO SCAN
## --------------------------
## by default, geantino (massless, chargeless virtual particle)
## events are generated uniformly of theta, phi
if  [ "$kopt" == 1 ]; then
	## output geant4_material_tracks.root
	## The result of the geantino scan will be a rootfile containing 
	## material tracks. Those contain the direction and production 
	## vertex of the geantino, the total material accumulated and all
	## the interaction points in the detector.
	${ACTS_BUILD_DIR}/bin/ActsExampleMaterialRecordingDD4hep  -j1 --dd4hep-input ${DETECTOR_XML_PATH}/${DETECTOR_NAME}.xml --output-root  -n ${nev} --gen-nparticles ${nparticles} #--gen-eta 3.5 4 --gen-phi -h
fi

## --------------------------
## ACTS mapping configuration
## --------------------------
if  [ "$kopt" == 2 ]; then
	# map ACTS gemoery to geometry-map.json
	${ACTS_BUILD_DIR}/bin/ActsExampleGeometryDD4hep -n1 -j1 --mat-output-file geometry-map  --dd4hep-input ${DETECTOR_XML_PATH}/${DETECTOR_NAME}.xml --output-json --mat-output-allmaterial true
	# convert geometry-map.json to config-map.json (to simplify modification)
	python3 ${ACTS_PATH}/Examples/Scripts/MaterialMapping/writeMapConfig.py geometry-map.json config-map.json
	echo ""
	echo "--------------DONE----------------"
	echo "config-map.json created. Please run materialmap_config.py next, or modify config-map.json manually."

	## Automatically edit the config-map.json to turn on mapMaterial on all approach 1 (right before representing surface)
	## ,approach 2 (right after representing surface) for endcaps
	## ,but NOT the central beampipe surface
	## Tested for Craterlake, many need to modify manually, and adjust binnings for other configurations
        # python materialmap_config.py
fi

## --------------------------
## Create map
## --------------------------

if  [ "$kopt" == 3 ]; then
	# turn config-map.json into modified geometry-map.json
	python3 ${ACTS_PATH}/Examples/Scripts/MaterialMapping/configureMap.py geometry-map.json config-map_new.json

	#----MAPPING------------
	# input: geant4_material_tracks.root, geometry-map.json
	# output: material-maps.json, root, cbor, materia-maps_tracks.root(from geantino)
	# ONLY map to surfaces not volumes for now
	${ACTS_BUILD_DIR}/bin/ActsExampleMaterialMappingDD4hep -j1 --input-root true --input-files geant4_material_tracks.root --mat-input-type file --mat-input-file geometry-map.json --output-root --output-json --output-cbor --mat-output-file material-maps --mat-mapping-surfaces true --mat-mapping-volumes false  --mat-mapping-read-surfaces false --dd4hep-input ${DETECTOR_XML_PATH}/${DETECTOR_NAME}.xml

	    # - ActsExampleMaterialMappingDD4hep -j 1 -n ${JUGGLER_N_EVENTS} --dd4hep-input ${DETECTOR_XML_PATH}/${DETECTOR_CONFIG}.xml --input-root true --input-files geant4_material_tracks.root --mat-input-type file --mat-input-file geometry-map.json --output-root --output-json --output-cbor --mat-output-file material-maps --mat-mapping-surfaces true --mat-mapping-read-surfaces false --mat-mapping-volumes true --mat-mapping-volume-stepsize 1


	#----Validation--------
	# output propagation-material.root
	${ACTS_BUILD_DIR}/bin/ActsExampleMaterialValidationDD4hep -n ${nev} --prop-ntests ${nparticles}  --mat-input-type file --mat-input-file material-maps.json --output-root --mat-output-file val-mat-map --dd4hep-input ${DETECTOR_XML_PATH}/${DETECTOR_NAME}.xml --prop-z0-sigma 0.0 --prop-d0-sigma 0.0
	# PropagationA   ERROR     Propagation reached the step count limit of 1000 (did 1000 steps)
	#  PropagationA   ERROR     Step failed with EigenStepperError:3: Step size adjustment exceeds maximum trials
	rm -rf Validation
	mkdir Validation

	root -l -b -q mat_map_local.C'("propagation-material.root","material-maps_tracks.root","Validation")'	# 

fi

## --------------------------
## More diagnostic plots
## --------------------------
# compare geantino scan v.s. material map surface by surface
if  [ "$kopt" == 4 ]; then

rm -rf Surfaces
mkdir Surfaces
cd Surfaces
mkdir prop_plot
mkdir map_plot
mkdir ratio_plot
mkdir dist_plot
mkdir 1D_plot
cd ..

root -l -b -q ${ACTS_PATH}/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_ratio.C'("propagation-material.root","material-maps_tracks.root",5000,"Surfaces/ratio_plot","Surfaces/prop_plot","Surfaces/map_plot")'
root -l -b -q ${ACTS_PATH}/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_dist.C'("material-maps_tracks.root",5000,"Surfaces/dist_plot")'
root -l -b -q ${ACTS_PATH}/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_1D.C'("material-maps_tracks.root",5000,"Surfaces/1D_plot")'

fi

