import os
import sys
import time
from liblibra_core import *
from libra_py.workflows.nbra import step3

# For excited states - Computing the excited states SDs and their overlaps and NACs
params_sd = {
          'lowest_orbital': 8192-100, 'highest_orbital': 8192+1000,
          'num_occ_states': 10, 'num_unocc_states': 10,
          'isUKS': 0, 'number_of_states': 10, 'tolerance': 0.01, 'verbosity': 0, 
          'use_multiprocessing': True, 'nprocs': 12, 
          'is_many_body': False, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': os.getcwd()+'/../step2/res',
          'logfile_directory': os.getcwd()+'/../step2/all_logfiles',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-mixed-basis-raw-ordered-identity',
          'outdir': os.getcwd()+'/res-mixed-basis-raw-ordered',
          'start_time': 900, 'finish_time': 3099, 'sorting_type': 'identity',
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }

step3.run_step3_sd_nacs_libint(params_sd)

