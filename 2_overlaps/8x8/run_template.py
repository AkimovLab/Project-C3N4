import os
import sys
from libra_py import CP2K_methods
from libra_py.workflows.nbra import step2


path = os.getcwd()
params = {}
# number of processors
params['nprocs'] = 25
# The CP2K and mpi executable
### Specific to CCR
params['mpi_executable'] = 'mpirun'
params['cp2k_exe'] = '/projects/academic/cyberwksp21/Software/cp2k-gnu/cp2k-v22.1/exe/local/cp2k.psmp'
params['istep'] = 
params['fstep'] = 
# Lowest and highest orbital, Here HOMO is 512
params['lowest_orbital'] = 8192-100
params['highest_orbital'] = 8192+1000
# extended tight-binding calculation type
params['isxTB'] = True
# DFT calculation type
params['isUKS'] = False
# Periodic calculations flag
params['is_periodic'] = True
# Set the cell parameters for periodic calculations
if params['is_periodic']:
    params['A_cell_vector'] = [113.9320831299   ,      0.0000000000    ,     0.0000000000]
    params['B_cell_vector'] = [  0.0017212471   ,     98.6744003146    ,     0.0000000000]
    params['C_cell_vector'] = [  0.0018975023   ,      0.0028002808    ,    14.9999996186]
    params['periodicity_type'] = 'XY'
    # Set the origin
    origin = [0,0,0]
    tr_vecs = params['translational_vectors'] = CP2K_methods.generate_translational_vectors(origin, [1,1,1],
                                                                                            params['periodicity_type'])
    
    print('The translational vectors for the current periodic system are:\n')
    print(tr_vecs)
    print(F'Will compute the S^AO between R(0,0,0) and {tr_vecs.shape[0]+1} translational vectors')

# The AO overlaps in spherical or Cartesian coordinates
params['is_spherical'] =  True
# Remove the molden files, which are large files for some systems, 
# after the computaion is done for tha system
params['remove_molden'] = True
# Cube visualization using VMD
params['cube_visualization'] = True
if params['cube_visualization']:
    # The only parts that we will change in this template are loading the cubes and rendering the images.
    params['vmd_input_template'] = '../vmd.tcl'
    params['states_to_plot'] = list(range(8183,8203))
    params['plot_phase_corrected'] = False
    params['vmd_exe'] = 'vmd'
    params['tachyon_exe'] = '/util/academic/vmd/1.9.2/lib/vmd/tachyon_LINUXAMD64'
    params['x_pixels'] = 3500
    params['y_pixels'] = 3500
    params['together_mode'] = True
    params['all_images'] = path + '/../all_images'
    params['image_format'] = 'bmp'
    params['remove_cube'] = True

# The results are stored in this folder
params['res_dir'] = path + '/../res'
params['all_pdosfiles'] = path + '/../all_pdosfiles'
params['all_logfiles'] = path + '/../all_logfiles'
# If the xTB calculations are needed, we need an OT procedure 
params['cp2k_ot_input_template'] = path + '/../es_ot_temp.inp'
params['cp2k_diag_input_template'] = path + '/../es_diag_temp.inp'
# The trajectory xyz file path
params['trajectory_xyz_filename'] = path + '/../../../1_molecular_dynamics/8x8/c3n4_8x8_MD_xTB-pos-1.xyz'

step2.run_cp2k_libint_step2(params)


