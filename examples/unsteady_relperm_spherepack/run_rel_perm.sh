#!bin/bash

# cp input/spheres4Palabos_105_100_100.dat tmp/

# cd /tmp

# matlab -r "addpath('/workspace/gigliotti/MultiphasePorousMediaPalabos_fork/MultiphasePorousMediaPalabos/post-processing/'); read_save_fluids; exit;"

# cd ../

mpirun -mca btl tcp,self -np 30 /workspace/gigliotti/MultiphasePorousMediaPalabos_fork/MultiphasePorousMediaPalabos/src/1-phase_LBM/permeability input_rel_perm.xml

mail -s 'Permeability stuff is really done!' alex.gigliotti@utexas.edu <<< 'Whoo! Rel Perms! :)'
