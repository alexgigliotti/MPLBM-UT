#!bin/bash

mpirun -mca btl tcp,self ../MultiphasePorousMediaPalabos/src/2-phase_LBM/ShanChen 2_phase_input.xml
mail -s 'Simulation is done!!' alex.gigliotti@utexas.edu <<< 'Yay! The simulation is done! :)'

