import numpy as np
import matplotlib.pyplot as plt
import h5py


# # Plotting n-states model timescales
plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
schemes = ['FSSH','IDA','mSDM']
colors = ['blue', 'orange', 'green']
lw = 1.
for c, scheme in enumerate(schemes):
    x = np.loadtxt(F'{scheme}.txt')
    if scheme=='IDA':
        label = 'ID-A'
    else:
        label = scheme
    line1 = plt.errorbar(range(7), x[:,0]/1000000, yerr=x[:,1]/1000000,
                         capsize=5, label=label,marker='s', lw=lw, ls='-', color=colors[c])
    print(F'Scheme {scheme}')
    
plt.ylabel('Timescale, ns')
plt.xlabel('C$_3$N$_4$ cell size')
# plt.yscale('log')
plt.xticks([0,1,2,3,4,5,6],['2x2','3x3','4x4','5x5','6x6','8x8','10x10'], fontsize=9)
plt.title('n-states')
plt.legend(fontsize=7.5)
plt.tight_layout()
plt.savefig(F'Fig3_id_n-states.jpg', dpi=600)


# # Plotting 2-states model timescales
plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
schemes = ['FSSH','IDA','mSDM']
colors = ['blue', 'orange', 'green']
lw = 1.
for c, scheme in enumerate(schemes):
    x = np.loadtxt(F'{scheme}-2.txt')
    if scheme=='IDA':
        label = 'ID-A'
    else:
        label = scheme
    line1 = plt.errorbar(range(7), x[:,0]/1000000, yerr=x[:,1]/1000000,
                         capsize=5, label=label,marker='s', lw=lw, ls='-', color=colors[c])
    print(F'Scheme {scheme}')
    
plt.ylabel('Timescale, ns')
plt.xlabel('C$_3$N$_4$ cell size')
# plt.yscale('log')
plt.xticks([0,1,2,3,4,5,6],['2x2','3x3','4x4','5x5','6x6','8x8','10x10'], fontsize=9)
plt.title('2-states')
plt.legend(fontsize=7.5)
plt.tight_layout()
plt.savefig(F'Fig3_id_2-states.jpg', dpi=600)


# # Plotting n-states and 2-states for each decoherence scheme
schemes = ['FSSH','IDA','mSDM']
colors = ['blue', 'orange', 'green']
lw = 1.
for c, scheme in enumerate(schemes):
    plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
    x = np.loadtxt(F'{scheme}.txt')
    line1 = plt.errorbar(range(7), x[:,0]/1000000, yerr=x[:,1]/1000000,
                         capsize=5, label='n-states',marker='s', lw=lw, ls='-')
    x = np.loadtxt(F'{scheme}-2.txt')
    line1 = plt.errorbar(range(7), x[:,0]/1000000, yerr=x[:,1]/1000000,
                         capsize=5, label='2-states',marker='o', lw=lw, ls='--')
    print(F'Scheme {scheme}')
    
    plt.ylabel('Timescale, ns')
    plt.xlabel('C$_3$N$_4$ cell size')
    # plt.yscale('log')
    plt.xticks([0,1,2,3,4,5,6],['2x2','3x3','4x4','5x5','6x6','8x8','10x10'], fontsize=9)
    if scheme=='IDA':
        title = 'ID-A'
    else:
        title = scheme
    plt.title(F'{title}')
    plt.legend(fontsize=7.5)
    plt.tight_layout()
    plt.savefig(F'Fig3_id_comp_active_space_{scheme}.jpg', dpi=600)

