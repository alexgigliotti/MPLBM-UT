#!bin/bash

mpirun -mca btl tcp,self ../MultiphasePorousMediaPalabos/src/1-phase_LBM/permeability rel_perm_input.xml
mail -s 'Permeability stuff is really done!' alex.gigliotti@utexas.edu <<< 'Whoo! Rel Perms! :)'
