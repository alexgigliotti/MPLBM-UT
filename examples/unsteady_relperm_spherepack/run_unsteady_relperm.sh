#! bin/bash

# Get image from Digital Rocks
bash fetch_image.sh

# Create 2-phase geometry
matlab -r 'create_geom_4_2phase; exit;'

# Run Shan-Chen
mpiexec -np 40 ../../src/2-phase_LBM/ShanChen input_spherepack.xml
mail -s 'Shan-Chen is done!!' alex.gigliotti@utexas.edu <<< 'Yay! The Shan-Chen is done! :)'

# Create rel perm geometries
matlab -r 'create_geoms_4_kr; exit;'

# Run rel perm simulation
# input folder = input/
# output folder = 4relperm/
mkdir 4relperm
mpiexec -np 40 ../../src/1-phase_LBM/permeability input_rel_perm.xml
mail -s 'Permeability stuff is really done!' alex.gigliotti@utexas.edu <<< 'Whoo! Rel Perms! :)'
