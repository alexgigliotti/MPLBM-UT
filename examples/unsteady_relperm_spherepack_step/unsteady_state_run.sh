###bash fetch_image.sh

# Create 2 phase geometry for flow simulation
#matlab -nodesktop -nosplash -r 'create_geom_4_2phase; exit;'

# Run Shan Chen LBM simulation
mpirun -np 8 ../../src/2-phase_LBM/ShanChen input_spherepack.xml

# Create rel perm geometry
#matlab -nodesktop -nosplash -r 'create_geoms_4_kr; exit;'

# Run rel perm simulation
##mkdir 4relperm
#mpirun -np 8 ../../src/1-phase_LBM/permeability input_rel_perm.xml
