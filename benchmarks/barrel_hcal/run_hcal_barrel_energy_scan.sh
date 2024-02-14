#!/bin/bash

export PARTICLE=$1
shift

E_file="sim_output/hcal_barrel_energy_scan_points_${PARTICLE}_${CI_JOB_ID}.txt"

if [ $# -gt 0 ] ; then
  E_VALUES=("$@")
else
  E_VALUES=(0.25 0.5 1 2 3 4 7 15 20)
fi

for E in ${E_VALUES[@]}
do
   export E_START="$E"
   export E_END="$E"
   export JUGGLER_FILE_NAME_TAG=hcal_barrel_${PARTICLE}_${E}
   bash benchmarks/barrel_hcal/run_hcal_barrel_particles.sh "${PARTICLE}" && echo "$E" >> "$E_file" || exit 1
   path_rootfiles="sim_output/energy_scan/${E}/"
   path_plots="results/energy_scan/${E}/"
   mkdir -p "$path_rootfiles"
   mkdir -p "$path_plots"
   ls -lthaR sim_output/
   mv "sim_output/sim_${JUGGLER_FILE_NAME_TAG}.edm4hep.root" "$path_rootfiles/sim_hcal_barrel_${PARTICLE}.edm4hep.root"
done

ls -lthaR sim_output
cat "$E_file"

exit 0
