#!/bin/bash

export PARTICLE=$1
E_file="sim_output/hcal_barrel_energy_scan_points_${PARTICLE}.txt"


#for E in 1 2 6 10
for E in 0.25 0.5 1 2 3 4 7 15 20
do
   export E_START="$E"
   export E_END="$E"
   bash benchmarks/barrel_hcal/run_hcal_barrel_particles.sh "${PARTICLE}" && echo "$E" >> "$E_file" || exit 1
   path_rootfiles="sim_output/energy_scan/${E}/"
   path_plots="results/energy_scan/${E}/"
   mkdir -p "$path_rootfiles"
   mkdir -p "$path_plots"
   ls -lthaR sim_output/
   mv "sim_output/sim_hcal_barrel_${PARTICLE}.root" "$path_rootfiles"
done

ls -lthaR sim_output
cat "$E_file"

exit 0
