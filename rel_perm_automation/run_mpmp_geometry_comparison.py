#!/workspace/gigliotti/auto_palabos_env/bin/python3.6

import numpy as np
import matplotlib.pyplot as plt
import os
import fileinput
import sys
from datetime import datetime

# Todo
# 1) Get this working
# 2) Write readme
# 3) Have everything in this file. Move bash file commands to this script
# 4) Create yaml input file for easy changes


class AutomatedMPMP():

    def __init__(self):

        # Get list of file names and links: geometry_links.txt
        links_and_names = 'geometry_names_and_links.txt'
        names_and_links_data = np.genfromtxt(links_and_names, delimiter="\t", dtype=str)
        self.geometry_names = names_and_links_data[:][:, 0]
        self.geometry_links = names_and_links_data[:][:, 1]

        # Set output directory
        self.output_directory = '/storage/gigliotti/geometry_comparison_test'

        # Size of the downloaded raw geometry (voxels per side)
        self.raw_geometry_size = 500

        # Fluid properties
        self.interfacial_tension = 36.23  # dynes/cm
        contact_angle = 0  # degrees
        self.contact_angle = contact_angle * np.pi/180  # Convert to radians

        # Set simulation size
        nx = 250 + 5  # Add 5 for padding on both ends of simulation (?)
        ny = 250
        nz = 250
        self.simulation_size = [nx, ny, nz]

        return

    def replace_line_in_file(self, file_to_edit, line_to_find_and_replace, replacement_line):

        search_and_replace_command = 'sed -i "/^' + line_to_find_and_replace + r"/c\\" + replacement_line + '" ' + file_to_edit
        os.system(search_and_replace_command)

        return
    
    def create_output_directory(self, i):
        geom_names = self.geometry_names
        output_dir = self.output_directory
        
        os.system("mkdir " + output_dir)  # create parent output directory
        os.system("mkdir " + output_dir + "/" + geom_names[i])  # create output for geometry
        os.system("mkdir " + output_dir + "/tmp")  # create output for simulation files if not there
        
        os.system("mkdir " + output_dir + "/tmp/2_phase_shanchen_output")
        os.system("mkdir " + output_dir + "/tmp/2_phase_shanchen_output/fluid_geometry_files")
        os.system("mkdir " + output_dir + "/tmp/2_phase_shanchen_output/fluid_geometry_files/input")
        os.system("mkdir " + output_dir + "/tmp/2_phase_shanchen_output/gifs")
        os.system("mkdir " + output_dir + "/tmp/2_phase_shanchen_output/numerical_outputs")

        os.system("mkdir " + output_dir + "/tmp/raw_geometry")

        os.system("mkdir " + output_dir + "/tmp/rel_perm_calculation")
        os.system("mkdir " + output_dir + "/tmp/rel_perm_calculation/fluid_geometry_input")
        os.system("mkdir " + output_dir + "/tmp/rel_perm_calculation/fluid_geometry_input/input")
        os.system("mkdir " + output_dir + "/tmp/rel_perm_calculation/output")
        os.system("mkdir " + output_dir + "/tmp/rel_perm_calculation/numerical_outputs")
        
        return
    
    def download_geometry_file(self, i):
        geom_names = self.geometry_names
        geom_links = self.geometry_links
        output_dir = self.output_directory
        
        os.system("wget -O " + output_dir + "/tmp/raw_geometry/" + geom_names[i] + ".raw " + geom_links[i])
        
        return 
    
    def update_and_run_geometry_creation_matlab_file(self, i):
        geom_names = self.geometry_names
        Nx = self.simulation_size[0]
        output_dir = self.output_directory
        geom_size = self.raw_geometry_size
       
        file_to_edit = "create_geom_4sim.m"
        line_to_find_and_replace = 'raw_geometry_file_to_open ='
        replacement_line = "raw_geometry_file_to_open = '" + output_dir + "/tmp/raw_geometry/" + geom_names[i] + ".raw';"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        line_to_find_and_replace = 'output_directory_for_porosity ='
        replacement_line = "output_directory_for_porosity = '" + output_dir + "/tmp/2_phase_shanchen_output/numerical_outputs/'"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        line_to_find_and_replace = 'geometry_file_output_name ='
        replacement_line = "geometry_file_output_name = '" + geom_names[i] + "_geometry';"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        line_to_find_and_replace = 'print_size ='
        replacement_line = "print_size = " + str(Nx-5) + "; %size of the Finneypack subset (in voxels per side)"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        line_to_find_and_replace = 'd_size ='
        replacement_line = "d_size = " + str(geom_size) + "; %voxels each side"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        line_to_find_and_replace = 'output_directory_for_geometry ='
        replacement_line = "output_directory_for_geometry = '" + output_dir + "/tmp/2_phase_shanchen_output/fluid_geometry_files/input/'"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        os.system('matlab -r "create_geom_4sim; exit"')
    
        return
    
    def update_2_phase_shan_chen_input_file(self, i):
        geom_names = self.geometry_names
        sim_size = self.simulation_size
        output_dir = self.output_directory

        # Update 2_phase_input.xml
        file_to_edit = '2_phase_input.xml'

        # Update geometry file location
        line_to_find_and_replace = '    <file_geom>'
        replacement_line = "    <file_geom> "+ output_dir + "/tmp/2_phase_shanchen_output/fluid_geometry_files/input/" + geom_names[i] + "_geometry_" + \
                           str(sim_size[0]) + "_" + str(sim_size[1]) + "_" + str(sim_size[2]) \
                           + ".dat" + " </file_geom>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        # Update simulation size
        line_to_find_and_replace = '    <size> <x>'
        replacement_line = "    <size> <x> " + str(sim_size[0]) + " </x> <y> " + str(sim_size[1]) + \
                           " </y> <z> " + str(sim_size[2]) + " </z> </size>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)
        
        # Update fluid 1 and 2 initial position
        line_to_find_and_replace = '    <fluid1> <x1>'
        replacement_line = "    <fluid1> <x1> 001 </x1> <y1> 001 </y1> <z1> 001 </z1> " + \
                           "<x2> 002 </x2> <y2> " + str(sim_size[1]) + " </y2> <z2> " + str(sim_size[2]) + " </z2> </fluid1>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        line_to_find_and_replace = '    <fluid2> <x1>'
        replacement_line = "    <fluid2> <x1> 003 </x1> <y1> 001 </y1> <z1> 001 </z1> " + \
                           "<x2> " + str(sim_size[0]) + " </x2> <y2> " + str(sim_size[1]) + " </y2> <z2> " + str(sim_size[2]) + " </z2> </fluid2>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        # Update output directory
        line_to_find_and_replace = '    <out_folder>'
        replacement_line = "    <out_folder> " + output_dir + "/tmp/2_phase_shanchen_output/ </out_folder>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        return
 
    def run_shan_chen_lbm_simulation(self):
        output_dir = self.output_directory
        os.system('bash 2_phase_auto.sh')

        # Organize outputs
        os.system('mv ' + output_dir + '/tmp/2_phase_shanchen_output/runnum.dat ' + output_dir + '/tmp/2_phase_shanchen_output/numerical_outputs/')
        os.system('mv ' + output_dir + '/tmp/2_phase_shanchen_output/rho_f1_*.gif ' + output_dir + '/tmp/2_phase_shanchen_output/gifs/')
        os.system('mv ' + output_dir + '/tmp/2_phase_shanchen_output/rho_f1_*.dat ' + output_dir + '/tmp/2_phase_shanchen_output/fluid_geometry_files/')
        os.system('mv ' + output_dir + '/tmp/2_phase_shanchen_output/porousMedium.vti ' + output_dir + '/tmp/2_phase_shanchen_output/fluid_geometry_files/')

        return

    def update_rel_perm_input_file(self, i):

        geom_names = self.geometry_names
        sim_size = self.simulation_size
        output_dir = self.output_directory

        # Update file name
        file_to_edit = 'rel_perm_input.xml'
        line_to_find_and_replace = '    <file_geom>'
        replacement_line = "    <file_geom> input/" \
                           + geom_names[i] + "_geometry_" + str(sim_size[0]) + "_" + str(sim_size[1]) + "_" + str(sim_size[2]) \
                           + " </file_geom>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        # Update simulation size
        file_to_edit = 'rel_perm_input.xml'
        line_to_find_and_replace = '    <size> <x>'
        replacement_line = "    <size> <x> " + str(sim_size[0]) + " </x> <y> " + str(sim_size[1]) + \
                           " </y> <z> " + str(sim_size[2]) + " </z> </size>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        # Update output folder
        file_to_edit = 'rel_perm_input.xml'
        line_to_find_and_replace = '    <out_f>'
        replacement_line = "    <out_f> " + output_dir + "/tmp/rel_perm_calculation/output/ </out_f>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        # Update input folder
        file_to_edit = 'rel_perm_input.xml'
        line_to_find_and_replace = '    <in_f>'
        replacement_line = "    <in_f> " + output_dir + "/tmp/rel_perm_calculation/fluid_geometry_input/ </in_f>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)
 
        # Get number of runs from Shan Chen Output. Multiply by 2 for wetting and non-wetting phases
        number_of_runs = np.loadtxt(output_dir + '/tmp/2_phase_shanchen_output/numerical_outputs/runnum.dat')
        self.number_of_runs = number_of_runs
        file_to_edit = 'rel_perm_input.xml'
        line_to_find_and_replace = '    <num>'
        replacement_line = "    <num> " + str(int(number_of_runs*2)) + " </num>"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        # Copy geometry file to fluid gemoetry directory per rel perm c++ file requirements
        os.system("cp " + output_dir + "/tmp/2_phase_shanchen_output/fluid_geometry_files/input/*.dat " \
                  + output_dir + "/tmp/rel_perm_calculation/fluid_geometry_input/input")

        return

    def create_geometry_for_rel_perms(self):
        output_dir = self.output_directory

        # Turn off percolation path (was causing problems before...)
        # medial_axis=0; %To find fluid medial axes and first percolation path
        file_to_edit = "../post-processing/read_save_fluids.m"
        line_to_find_and_replace = 'medial_axis='
        replacement_line = "medial_axis=0; %To find fluid medial axes and first percolation path"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        # Update directory
        # tmp/2_phase_shanchen_output/       directory=
        file_to_edit = "../post-processing/read_save_fluids.m"
        line_to_find_and_replace = 'directory='
        replacement_line = "directory= '" + output_dir + "/tmp/2_phase_shanchen_output/fluid_geometry_files/';"
        self.replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line)

        # Generate geometry files for rel perm calculation
        add_path = "addpath('../post-processing');"
        os.system('matlab -r " diary error.txt;' + add_path +' read_save_fluids; exit"')

        # Organize output
        os.system('mv ' + output_dir + '/tmp/2_phase_shanchen_output/fluid_geometry_files/*.csv ' + output_dir + '/tmp/2_phase_shanchen_output/fluid_geometry_files/Sw.csv')
        os.system('mv ' + output_dir + '/tmp/2_phase_shanchen_output/fluid_geometry_files/lattice_f* ' + output_dir + '/tmp/rel_perm_calculation/fluid_geometry_input/')
        os.system('mv ' + output_dir + '/tmp/2_phase_shanchen_output/fluid_geometry_files/Sw.csv ' + output_dir + '/tmp/2_phase_shanchen_output/numerical_outputs/')

        return

    def run_rel_perm_simulation(self):

        output_dir = self.output_directory
        os.system('bash rel_perm_auto.sh')
        os.system('mv perm_data.csv relperm_data.csv ' + output_dir + '/tmp/rel_perm_calculation/numerical_outputs/')

        return

    def save_results_and_clear_temp(self, i):

        geom_names = self.geometry_names
        output_dir = self.output_directory

        # Copy all results to the output directory for the geometry being simulated
        os.system("cp -a " + output_dir + "/tmp/. " + output_dir + "/" + geom_names[i] + "/")

        # Delete temporary files and geometry input file in automation directory
        os.system("find " + output_dir + "/tmp/ -type f -delete")
        # os.system("find input/ -type f -delete")

        return
 
    def parse_simulation_data(self, i):

        geom_names = self.geometry_names
        output_dir = self.output_directory

        shan_chen_outputs_dir = output_dir + "/" + geom_names[i] + "/2_phase_shanchen_output/numerical_outputs/"
        rel_perm_outputs_dir = output_dir + "/" + geom_names[i] + "/rel_perm_calculation/numerical_outputs/"

        number_of_runs = int(np.loadtxt(shan_chen_outputs_dir + 'runnum.dat'))
        self.number_of_runs = number_of_runs
        num_runs = number_of_runs

        self.water_saturation = np.genfromtxt(shan_chen_outputs_dir + 'Sw.csv', delimiter=',')
        self.porosity = np.genfromtxt(shan_chen_outputs_dir + 'porosity.csv', delimiter=',')
        # self.capillary_number = np.genfromtxt(shan_chen_outputs_dir + 'capillary_number_data.csv', delimiter=',')[0:num_runs]
        # self.mobility_ratio_data = np.genfromtxt(shan_chen_outputs_dir + 'mobility_ratio_data.csv', delimiter=',')

        # self.psudo_capillary_pressure = np.genfromtxt(shan_chen_outputs_dir + 'capillary_pressure_data.csv', delimiter=',')[0:num_runs]
        Gc = 0.9
        rho_w = 2  # density wetting phase
        rho_nw = np.linspace(1, 2, num_runs)  # density nonwetting phase
        self.capillary_pressure = (rho_w + rho_nw)/3 + (Gc*rho_w*rho_nw/3)

        self.absolute_permeability = np.genfromtxt(rel_perm_outputs_dir + 'perm_data.csv', delimiter=',')[0]
        self.relative_permeability = np.array(np.genfromtxt(rel_perm_outputs_dir + 'relperm_data.csv', delimiter=','))
        self.relative_permeability_wetting = self.relative_permeability[0:num_runs]
        self.relative_permeability_nonwetting = self.relative_permeability[num_runs:-1]

        return

    def plot_results(self, i):

        geom_names = self.geometry_names
        output_dir = self.output_directory
        num_runs = self.number_of_runs

        Sw = self.water_saturation
        phi = self.porosity
        sigma = self.interfacial_tension
        theta = self.contact_angle
        # Ca = self.capillary_number
        # M = self.mobility_ratio_data
        Pc = self.capillary_pressure
        k = self.absolute_permeability
        krw = self.relative_permeability_wetting
        krnw = self.relative_permeability_nonwetting

        # Create output directory for figures
        figures_dir = output_dir + "/" + geom_names[i] + "/figures"
        os.system("mkdir " + figures_dir)
        
        # Replace underscores in names with a space so it looks nice
        geom_names[i].replace("_", " ")

        # Normalize capillary pressure by Leverett J function
        J = Pc*np.sqrt(k/phi) / (sigma * np.abs(np.cos(theta)))

        # Capillary pressure curve
        plt.figure()
        plt.plot(Sw, J, '-o')
        plt.xlabel("Wetting Phase Saturation")
        plt.ylabel("Capillary Pressure, J-Function")
        plt.title(geom_names[i] + " LBM Drainage Capillary Pressure Curve, J-Function")
        plt.grid()
        plt.savefig(figures_dir + "/capillary_pressure_curve.png")

        # Rel perm curves
        plt.figure()
        plt.plot(Sw, krw, '-o')
        plt.plot(Sw, krnw, '-o')
        plt.xlabel("Wetting Phase Saturation")
        plt.ylabel("Relative Permeability")
        plt.title(geom_names[i] + " LBM Drainage Relative Permeability Curve")
        plt.legend(['krw', 'krnw'])
        plt.grid()
        plt.savefig(figures_dir + "/rel_perm_curve.png")

        return

    def email_when_done(self, i, end_time):
        output_dir = self.output_directory
        geom_names = self.geometry_names

        # mail -s 'Simulation is done!!' alex.gigliotti@utexas.edu <<< 'Yay! The simulation is done! :)'
        subject = 'Everything is done for ' + geom_names[i] + '!!'
        email = 'alex.gigliotti@utexas.edu'
        attachment_1 = output_dir + "/" + geom_names[i] + '/figures/capillary_pressure_curve.png'
        attachment_2 = output_dir + "/" + geom_names[i] + '/figures/rel_perm_curve.png'
        message = 'Yay! It took this amount of time to run: ' + str(end_time) + '\n\n :D'
        # print("mail -s " + "'" + subject + "' -a " + attachment_1 + ' -a ' + attachment_2 + ' ' + email + " <<< " + "'" + message + "'")
        os.system("mail -s " + "'" + subject + "' -a " + attachment_1 + ' -a ' + attachment_2 + ' ' + email + " <<< " + "'" + message + "'")

        return


start_time = datetime.now()
mpmp = AutomatedMPMP()

for i in range(len(mpmp.geometry_names)):

    # Start timer
    start_time = datetime.now()
    
    # 1) create output directory
    mpmp.create_output_directory(i)
    
    # 2) download and extract geometry file to tmp/input
    mpmp.download_geometry_file(i)
    
    # 3) Update file names in matlab file
    mpmp.update_and_run_geometry_creation_matlab_file(i)
    
    # 4) Edit xml input for Shan Chen
    mpmp.update_2_phase_shan_chen_input_file(i)

    # 5) Run Shan Chen
    mpmp.run_shan_chen_lbm_simulation()
 
    # 6) Update xml input file for rel perms
    mpmp.update_rel_perm_input_file(i)

    # 7) Generate geometries for rel perms
    mpmp.create_geometry_for_rel_perms()    

    # 8) Run Rel Perm simulation
    mpmp.run_rel_perm_simulation()

    # 9) Save results and clear tmp/ for next simulation
    mpmp.save_results_and_clear_temp(i)

    # 10) Parse data and create plots
    mpmp.parse_simulation_data(i)
    mpmp.plot_results(i)

    # 11) Email when done with the time!
    end_time = datetime.now() - start_time
    mpmp.email_when_done(i, end_time)

    break  # delete after testing!!

