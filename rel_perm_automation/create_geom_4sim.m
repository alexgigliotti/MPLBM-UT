% creates a geometry for simulation

%% opening the file
addpath ('../pre-processing') %pre-precesing libraries
d_size = 500; %voxels each side
raw_geometry_file_to_open = '/storage/gigliotti/geometry_comparison_test/tmp/raw_geometry/100_finney_pack_test.raw';
f1 = fopen(raw_geometry_file_to_open,'r'); %read raw file
fp = fread(f1, d_size*d_size*d_size,'uint8=>uint8');
fp = reshape(fp, d_size,d_size,d_size);

%% selecting a smaller subset
print_size = 100; %size of the Finneypack subset (in voxels per side)
fp_printing = fp(1:print_size, 1:print_size, 1:print_size);
figure();imagesc(fp_printing(:,:,uint8(print_size/2)));
title('Cross-section of the simulation subset')

%% eliminating non-connected regions 
connect = 6; % pixel connectivity 6, 18, 26
output_directory_for_porosity = '/storage/gigliotti/geometry_comparison_test/tmp/2_phase_shanchen_output/numerical_outputs/'
fp = eliminate_isolatedRegions(fp, connect, output_directory_for_porosity); %for better convergence

%% making a computationally efficent domain for sim
geometry_file_output_name = '100_finney_pack_test_geometry';
name = [geometry_file_output_name];
add_mesh   = true; % add a neutral-wet mesh at the end of the domain
num_slices = 2;    % add n empty slices at the beggining and end of domain 
                   % for pressure bcs
swapXZ = true;     % Swap x and z data if needed to ensure Palabos simulation in Z-direction
scale_2 = false;   % Double the grain (pore) size if needed to prevent single pixel throats
                   % for tight/ low porosity geometries
output_directory_for_geometry = '/storage/gigliotti/geometry_comparison_test/tmp/2_phase_shanchen_output/fluid_geometry_files/input/'

palabos_3Dmat = create_geom_edist(fp_printing,name,num_slices, add_mesh, ...
                                  swapXZ, scale_2, output_directory_for_geometry);
                                    %provides a very  efficient
                                    %geometry for Palabos


% old version
% palabos_3Dmat = mat2dat_4lbm(fp_printing,name,1); %although this function is slow, it still works

%% Mixed Wettability (the user could experiment with this)
rng(123)
rnd_array = rand(size(palabos_3Dmat) );

palabos_3Dmat_mixedWet = palabos_3Dmat;
palabos_3Dmat_mixedWet(palabos_3Dmat==1 & rnd_array>0.5)=3;

