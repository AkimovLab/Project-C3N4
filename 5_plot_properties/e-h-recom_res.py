import os
import sys
import glob
import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit
from libra_py import units
import matplotlib.pyplot as plt


def linear_func(t, tau):
    return 1-t/tau
def exp_func(t, tau):
    return np.exp(-t/tau)


def compute_z_val(n):
    if n==0:
        z = 0.000
    elif n==1:
        z = 12.706
    elif n==2:
        z = 4.303
    elif n==3:
        z = 3.182
    elif n==4:
        z = 2.776
    elif n==5:
        z = 2.571
    elif n==6:
        z = 2.447
    elif n==7:
        z = 2.365
    elif n==8:
        z = 2.306
    elif n==9:
        z = 2.262
    elif n==10:
        z = 2.228
    elif n==11:
        z = 2.201
    elif n==12:
        z = 2.179
    elif n==13:
        z = 2.160
    elif n==14:
        z = 2.145
    elif n==15:
        z = 2.131
    elif n==16:
        z = 2.120
    elif n==17:
        z = 2.110
    elif n==18:
        z = 2.101
    elif n==19:
        z = 2.093
    elif n==20:
        z = 2.086
    elif n==21:
        z = 2.080
    elif n==22:
        z = 2.074
    elif n==23:
        z = 2.069
    elif n==24:
        z = 2.064
    elif n==25:
        z = 2.060
    elif n==26:
        z = 2.056
    elif n==27:
        z = 2.052
    elif n==28:
        z = 2.048
    elif n==29:
        z = 2.045
    elif n==30:
        z = 2.042
    elif n>30 and n<=40:
        z = 2.021
    elif n>40 and n<=50:
        z = 2.009
    elif n>50 and n<=60:
        z = 2.000
    elif n>60 and n<=80:
        z = 1.990
    elif n>80 and n<=100:
        z = 1.984
    elif n>100 and n<=120:
        z = 1.980
    elif n>120:
        z = 1.960
    return z


# In[20]:


dt = 1.0
R2 = 0.1
sizes = ['2x2','3x3','4x4','5x5','6x6','8x8','10x10']
schemes = ['FSSH', 'IDA', 'mSDM'] #,'DISH']
for size in sizes:
    print(F'-----------------Fitting data for size {size} -----------------------')
    for scheme in schemes:
        print(F'----------------{scheme}------------')
        if scheme=='IDA':
            title = 'ID-A'
        else:
            title = scheme
        md_time = np.arange(0,10000,dt)
        # The power - we do not use that but I'll keep it 
        for i in [0]:
            taus = []
            plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
            files = glob.glob(F'../4_namd/{size}/0_namd_res_*/{scheme}/SH_pop*txt')  
            #print('Files for i ', i)
            c = 0
            c1 = 0
            min_val = 1.0
            for file in files:
                x = np.loadtxt(file)
                if x.shape[0]==10000:
                    c += 1
                    if c==1:
                        x_ave = np.zeros(x.shape)
                    popt, pcov = curve_fit(exp_func, md_time, 1-x[:,0],bounds=([0.0],[40000000]))
                    tau = popt[0]
                    residuals  = 1-x[:,0] - exp_func(md_time, *popt)
                    ss_res     = np.sum(residuals**2)
                    ss_tot     = np.sum((1-x[:,0] - np.mean(1-x[:,0]))**2)
                    r_squared  = 1.0 - (ss_res / ss_tot)
                    if r_squared>R2:
                        c1 += 1
                        taus.append(tau)
                        plt.plot(md_time, 1-x[:,0], color='gray')
                        print('R2:',r_squared,', tau:', tau/1000000,'ns')
                        x_ave += x
                        min_val = min([min_val, np.min(1-x[:,0])])
            taus = np.array(taus)
            N = len(taus)
            Z = compute_z_val(N)
            s = np.std(taus)
            tau_ave = np.average(taus)
            error_bar = Z*s/np.sqrt(N)
            print(F'Average timscale for power {i}:', tau_ave/1000000, ' +- ',error_bar/1000000, ' ns, number of fitted samples:', c1)


            print(F"Special flag for grepping the data! size {size} and {scheme} average time scale: {tau_ave} +- {error_bar} fs")


            if c1>0:
                plt.plot(exp_func(md_time, tau_ave), color='r')
                plt.plot(exp_func(md_time, tau_ave+error_bar), ls='--', color='r')
                plt.plot(exp_func(md_time, tau_ave-error_bar), ls='--', color='r')
            if scheme=='IDA':
                plt.title(F'{size}, ID-A')
            else:
                plt.title(F'{size}, {scheme}')
            if tau_ave > 1000000:
                plt.text(-20, min_val,'%0.2f$\pm$%0.2f ns'%(tau_ave/1000000, error_bar/1000000), fontsize=14, fontweight='bold')
            elif tau_ave > 1000.0 and tau_ave < 1000000:
                plt.text(-20, min_val,'%0.2f$\pm$%0.2f ps'%(tau_ave/1000, error_bar/1000), fontsize=14, fontweight='bold')
            elif tau_ave < 1000:
                plt.text(-20, min_val,'%0.2f$\pm$%0.2f fs'%(tau_ave, error_bar), fontsize=14, fontweight='bold')
            plt.xlabel('Time, fs')
            plt.ylabel('Population')
            plt.tight_layout()
            plt.savefig(F'{size}_{scheme}_{i}_R2_{R2}_type2_2states.jpg', dpi=600)
            plt.close()

