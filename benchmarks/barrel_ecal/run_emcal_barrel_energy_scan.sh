#!/bin/bash

export PARTICLE=$1
E_file="sim_output/emcal_barrel_energy_scan_points_${PARTICLE}.txt"

MIN_N_EVENTS=750
if (( $JUGGLER_N_EVENTS < $MIN_N_EVENTS )); then
   ORIGINAL_JUGGLER_N_EVENTS=$JUGGLER_N_EVENTS 
   export JUGGLER_N_EVENTS=$MIN_N_EVENTS
   echo "Setting JUGGLER_N_EVENTS to ${MIN_N_EVENTS}"
fi
#0.5 1 2 3 4 7 15 20
for E in 0.5 1 2 5 10
do
   export E_START="$E"
   export E_END="$E"   
   bash benchmarks/barrel_ecal/run_emcal_barrel_particles.sh "${PARTICLE}" && echo "$E" >> "$E_file" || exit 1
   path_rootfiles="sim_output/energy_scan/${E}/"
   path_plots="results/energy_scan/${E}/"
   mkdir -p "$path_rootfiles"
   mkdir -p "$path_plots"
   ls -lthaR sim_output/
   mv "sim_output/sim_emcal_barrel_${PARTICLE}.root" "$path_rootfiles"
done

export JUGGLER_N_EVENTS=$ORIGINAL_JUGGLER_N_EVENTS

ls -lthaR sim_output
cat "$E_file"

exit 0
