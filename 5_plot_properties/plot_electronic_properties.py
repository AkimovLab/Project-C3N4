#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import time
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from liblibra_core import *
from libra_py.workflows.nbra import step3
from libra_py import units, data_stat, data_io
import h5py


# # Reading the energies 

# In[2]:


# Generate the energies
sizes = ['2x2','3x3','4x4','5x5','6x6','8x8','10x10']
for size in sizes:
    energy_files = glob.glob(F'../3_nacs/{size}/res-mixed-basis-raw-ordered-identity/Hvib_sd*re*')
    energy_files = data_io.sort_hvib_file_names(energy_files)
    # print('Sorted energy files are:', energy_files)
    dt = 1.0 # fs
    energies = []
    #for file in energy_files:
    #    diag_elements = np.diag(sp.load_npz(file).todense().real)
    #    energies.append(diag_elements)
    #energies = np.array(energies)*units.au2ev
    
    # Lazy one
    for i in range(len(energy_files)):
        diag_elements = np.diag(sp.load_npz(energy_files[i]).todense().real)
        energies.append(diag_elements)
    energies = np.array(energies)*units.au2ev
    
    
    # # Plotting the energies vs time
    
    # In[3]:
    
    
    plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
    md_time = np.arange(0,energies.shape[0])
    print(energies.shape)
    for i in range(energies.shape[1]):
        plt.plot(md_time, energies[:,i]-energies[:,0],lw=1.0)
    plt.title(F'{size} C$_3$N$_4$, Mixed basis', fontsize=10)
    plt.ylabel('Energy, eV', fontsize=10)
    plt.xlabel('Time, fs', fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig(F'{size}_E_vs_Time_identity.jpg', dpi=600)
    
    
    # # Reading NACs, saving the distribution as HDF5 files, and computing the average NAC map for plot
    # 
    # 
    # The reason we need to save the distribution data as HDF5 is because we want to plot them all together in another script.
    
    # In[4]:
    
    
    # Here we want to find the index for which we have the average 1.0 eV above LUMO
    energies_ave = np.average(energies,axis=0)
    energies_ave -= energies_ave[1]
    print('The average energies are:\n',energies_ave)
    max_indx = np.min(np.where(energies_ave>0.2)[0])
    print('The index is:', max_indx)
    
    
    # In[ ]:
    
    
    # %matplotlib notebook
    # Plot NAC map
    # Choose this active space only up to index which has 1.0 eV above LUMO
    a = list(range(max_indx))
    # print(a)
    nac_files = glob.glob(F'../3_nacs/{size}/res-mixed-basis-raw-ordered/Hvib_sd*im*')
    # print(nac_files)
    nac = []
    for c, nac_file in enumerate(nac_files):
        # Choosing only the indices that are within the 1.0 eV excitation energy
        nac_mat = sp.load_npz(nac_file).todense().real[a,:][:,a]
        if c==0:
            nac_ave = np.zeros(nac_mat.shape)
        nac_ave += np.abs(nac_mat)
        # The distribution of the NACs
        for i in range(nac_mat.shape[0]):
            for j in range(nac_mat.shape[0]):
                if j != i:
                    nac_ij = np.abs(nac_mat[i,j])* 1000.0 * units.au2ev
                    x_mb = MATRIX(1,1)
                    x_mb.set(0, 0, nac_ij )
                    nac.append( x_mb )
    bin_supp, dens, cum = data_stat.cmat_distrib( nac, 0, 0, 0, 0, 500, 0.1)
    # plt.plot(bin_supp, dens)
    # Saving the data so that we can plot it with other supercells NAC distributions
    with h5py.File(F"{size}_nac_dist.hdf", "w") as f:
        #g = f.create_group("data")
        f.create_dataset("bin_supp", data = bin_supp)
        f.create_dataset("dens", data = dens)
    
    nac_ave *= 1000*units.au2ev/c
    nstates = nac_ave.shape[0]
    np.savetxt('nac_ave', nac_ave)
    #sys.exit(0)
    # # Testing the HDF5 file for plot
    
    # In[ ]:
    
    # # Plotting the average NAC map
    
    # In[ ]:
    
    
    plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
    plt.imshow(np.flipud(nac_ave), cmap='hot', extent=(0,nstates,0,nstates))#, vmin=0, vmax=150)
    plt.xlabel('State index', fontsize=10)
    plt.ylabel('State index', fontsize=10)
    plt.colorbar().ax.set_title('meV', fontsize=8)
    plt.title(F'{size} C$_3$N$_4$ - Mixed basis', fontsize=10)
    plt.tight_layout()
    # Some data may be needed in the manuscript
    for i in range(1,len(nac_ave[0])):
        print(F'The average NAC between the ground and the {i}th excited-state is',nac_ave[0,i], 'meV')
    plt.savefig(F'{size}_c3n4_mixed_basis_basis_nac_map.jpg', dpi=600)
    
    # # Plotting the PDOS
    
    # In[ ]:
    
    
    def gaussian_function(a, mu, sigma, num_points, x_min, x_max):
        pre_fact = (a/sigma)/(np.sqrt(2*np.pi))
        x = np.linspace(x_min, x_max, num_points)
        x_input = np.array((-1/2)/(np.square(sigma))*np.square(x-mu))
        gaussian_fun = pre_fact*np.exp(x_input)
        
        return x, gaussian_fun
        
    def gaussian_function_vector(a_vec, mu_vec, sigma, num_points, x_min, x_max):
        for i in range(len(a_vec)):
            if i==0:
                sum_vec = np.zeros(num_points)
            energy_grid, conv_vec = gaussian_function(a_vec[i], mu_vec[i], sigma, num_points, x_min, x_max)
            sum_vec += conv_vec
        return energy_grid, sum_vec
    
    
    # In[ ]:
    
    
    plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
    path_to_all_pdos = os.getcwd()+F'/../2_overlaps/{size}/all_pdosfiles'
    atoms = ['C', 'N']
    # Original
    # orbitals_cols = [[3], range(4,7), range(7,12), range(12,19)]
    # orbitals = ['s', 'p','d','f']
    
    # Added by 10/14/2022
    orbitals_cols = [[3], [4], [5], [6], range(7,12), range(12,19)]
    orbitals = ['s','py', 'pz', 'px','d','f']
    
    npoints = 3000
    sigma = 0.05
    shift = 3.0 # eV
    ave_pdos_convolved_all = []
    for c1,i in enumerate([1,2]):
        pdos_files = glob.glob(path_to_all_pdos+F'/*k{i}*.pdos')
        for c2, pdos_file in enumerate(pdos_files):
            pdos_mat = np.loadtxt(pdos_file)
            if c2==0:
                pdos_ave = np.zeros(pdos_mat.shape)
            pdos_ave += pdos_mat
        pdos_ave /= c2+1
        pdos_ave[:,1] *= units.au2ev
        e_min = np.min(pdos_ave[:,1])-shift
        e_max = np.max(pdos_ave[:,1])+shift
        homo_level = np.max(np.where(pdos_ave[:,2]==2.0))
        homo_energy = pdos_ave[:,1][homo_level]
        for c3, orbital_cols in enumerate(orbitals_cols):
            try:
                sum_pdos_ave = np.sum(pdos_ave[:,orbital_cols],axis=1)
                ave_energy_grid, ave_pdos_convolved = gaussian_function_vector(sum_pdos_ave, pdos_ave[:,1], sigma,
                                                                                   npoints, e_min, e_max)
                ave_pdos_convolved_all.append(ave_pdos_convolved)
                pdos_label = atoms[c1]+F', {orbitals[c3]}'
                plt.plot(ave_energy_grid-homo_energy, ave_pdos_convolved, label=pdos_label)
            except:
                pass
        #try:
        #    sum_pdos_ave = np.sum(pdos_ave[:,3::],axis=1)
        #    ave_energy_grid, ave_pdos_convolved = gaussian_function_vector(sum_pdos_ave, pdos_ave[:,1], sigma,
        #                                                                       npoints, e_min, e_max)
        #    ave_pdos_convolved_all.append(ave_pdos_convolved)
        #    pdos_label = atoms[c1]
        #    plt.plot(ave_energy_grid-homo_energy, ave_pdos_convolved, label=pdos_label)
        #except:
        #    pass
    
    
    ave_pdos_convolved_total = np.sum(np.array(ave_pdos_convolved_all),axis=0)
    plt.plot(ave_energy_grid-homo_energy, ave_pdos_convolved_total, color='black', label='Total')
    plt.legend(fontsize=7)
    plt.xlim(-5,5)
    plt.xticks([-5,-2.5,0.0,2.5,5], fontsize=10)
    plt.ylabel('pDOS, 1/eV', fontsize=10)
    plt.xlabel('Energy, eV', fontsize=10)
    plt.title(F'{size} C$_3$N$_4$, 300 K', fontsize=10)
    plt.tight_layout()
    # plt.savefig('../C3N4_pdos.jpg', dpi=600)
    plt.savefig(F'{size}_C3N4_pdos.jpg', dpi=600)



ax = plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
lw=2.0
for i in [2,3,4,5,6,8,10]:
    print(i)
    hf = h5py.File(F'{i}x{i}_nac_dist.hdf','r')
    bin_supp = hf.get('bin_supp')[:]
    dens = hf.get('dens')[:]
    plt.plot(bin_supp, dens, label=F'{i}x{i}',lw=lw)

# For 0<x<1
plt.yscale('log')
plt.xlim(0,1)
plt.xlabel('$\hbar$|NAC|, meV')
plt.ylabel('PD, 1/meV')
plt.tight_layout()
plt.savefig('nac_dist_1.jpg', dpi=600)
# For 0<x<100 but 0.0001<y<0.1
plt.xlim(0,100)
plt.ylim(1e-4,0.1)
plt.tight_layout()
plt.savefig('nac_dist_2.jpg', dpi=600)


