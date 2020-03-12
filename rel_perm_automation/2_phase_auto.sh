#!bin/bash

mpirun -mca btl tcp,self -np 10 ../src/2-phase_LBM/ShanChen 2_phase_input.xml
mail -s 'Shan-Chen is done!!' alex.gigliotti@utexas.edu <<< 'Yay! The Shan-Chen is done! :)'

