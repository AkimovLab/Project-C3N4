import os
import sys
from libra_py import CP2K_methods


run_slurm = True
submit_template = 'submit_template.slm'
run_python_file = 'run_template.py'
istep = 900
fstep = 3998
njobs = 20

# Removing the previous folders if existed. You can keep them as well 
# but Libra will overwrite some of the data if their names are the same
#os.system('rm -rf res job* all_logfiles all_pdosfiles')

print('Distributing jobs...')
CP2K_methods.distribute_cp2k_libint_jobs(submit_template, run_python_file, istep, fstep, njobs, run_slurm)

