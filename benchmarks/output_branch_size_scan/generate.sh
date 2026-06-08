#!/bin/bash
set -Euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR
IFS=$'\n\t'

NUM_EVENTS=400
INPUT_FILE=root://dtn-eic.jlab.org//work/eic2/EPIC/EVGEN/DIS/NC/18x275/minQ2=1/pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.hepmc3.tree.root

DETECTOR_CONFIG=epic_craterlake
EBEAM=18
PBEAM=275

npsim \
    --runType batch \
    --random.seed 1 \
    --random.enableEventSeed \
    --printLevel WARNING \
    --skipNEvents 0 \
    --numberOfEvents 400 \
    --filter.tracker 'edep0' \
    --hepmc3.useHepMC3 true \
    --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}${EBEAM:+${PBEAM:+_${EBEAM}x${PBEAM}}}.xml \
    --inputFiles ${INPUT_FILE} \
    --outputFile current_campaign.edm4hep.root
    
eicrecon \
    -Ppodio:output_file="current_campaign.eicrecon.tree.edm4eic.root" \
    -Pjana:warmup_timeout=0 -Pjana:timeout=0 \
    -Pplugins=janadot \
    "current_campaign.edm4hep.root"


