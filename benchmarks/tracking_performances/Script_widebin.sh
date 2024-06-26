#!/bin/bash
# Script; shyam kumar; INFN Bari, Italy
# shyam.kumar@ba.infn.it; shyam055119@gmail.com
rm -rf truthseed/ realseed/ *.root
mkdir -p truthseed/pi-/mom realseed/pi-/mom truthseed/pi-/dca realseed/pi-/dca
mom_array=(0.5 1.0 2.0 5.0 10.0 15.0)
particle_array=("pi-")
filename=("tracking_output") 
etabin_array=(-3.5 -2.5 -1.0 1.0 2.5 3.5)
nevents=10000

# run the simulation
source ../epic/install/bin/thisepic.sh 
dd_web_display --export -o epic_craterlake_tracking_only.root ../epic/install/share/epic/epic_craterlake_tracking_only.xml
for ((i=0; i<${#mom_array[@]}; i++)); do
npsim --compactFile ../epic/install/share/epic/epic_craterlake_tracking_only.xml --outputFile sim${mom_array[i]}.edm4hep.root --numberOfEvents $nevents --enableGun --gun.thetaMin 3*deg --gun.thetaMax 177*deg --gun.distribution eta --gun.particle pi- --gun.momentumMin ${mom_array[i]}*GeV --gun.momentumMax ${mom_array[i]}*GeV --gun.multiplicity 1 --random.seed 100000
done
# run the reconstruction
source ../epic/install/setup.sh 
source ../EICrecon/install/bin/eicrecon-this.sh

for ((i=0; i<${#mom_array[@]}; i++)); do
eicrecon \
 -Pnthreads=1 \
 -Pjana:debug_plugin_loading=1 \
 -Pjana:nevents=$nevents \
 -Pacts:MaterialMap=calibrations/materials-map.cbor \
 -Ppodio:output_file="${filename}"_${mom_array[i]}.edm4eic.root \
 -Pdd4hep:xml_files=../epic/install/share/epic/epic_craterlake_tracking_only.xml   \
 -Ppodio:output_collections="MCParticles,CentralCKFTrajectories,CentralCKFTrackParameters,CentralCKFSeededTrackParameters,CentralTrackVertices" \
 sim${mom_array[i]}.edm4hep.root
done 

# run the analysis
for ((iparticle=0; iparticle<${#particle_array[@]}; iparticle++)); do
# truth seeding
for ((i=0; i<${#mom_array[@]}; i++)); do
root -b -l -q Tracking_Performances.C'("'${filename}'","'${particle_array[iparticle]}'",'${mom_array[i]}',0.15,"")'
done

# real seeding
for ((i=0; i<${#mom_array[@]}; i++)); do
root -b -l -q Tracking_Performances.C'("'${filename}'","'${particle_array[iparticle]}'",'${mom_array[i]}',0.15,"Seeded")'
done
done
cd truthseed/pi-/dca
hadd final_hist_dca_truthseed.root *.root
cd ../../../realseed/pi-/dca/
hadd final_hist_dca_realseed.root *.root
cd ../../../

rm -rf Final_Results/ Debug_Plots/ 
mkdir -p Final_Results/pi-/mom Final_Results/pi-/dca  Debug_Plots/truth/pi-/mom  Debug_Plots/truth/pi-/dca  Debug_Plots/real/pi-/mom Debug_Plots/real/pi-/dca
# loop over particles
for ((iparticle=0; iparticle<${#particle_array[@]}; iparticle++)); do
for ((i=0; i<${#etabin_array[@]}-1; i++)); do
xmax_hist=0.3 
if [ $i == 2 || $i == 1 ]; then 
xmax_hist=0.01 
fi
root -b -l -q doCompare_truth_real_widebins_mom.C'("'${particle_array[iparticle]}'",'${etabin_array[i]}','${etabin_array[i+1]}','$xmax_hist')'
done
done

