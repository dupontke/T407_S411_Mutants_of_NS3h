#!/Users/dupontke/anaconda3/bin/python3
# ----------------------------------------
# USAGE:

# ./wt_helicase_time_data.py system comp_cells buffer

# ----------------------------------------
# PREAMBLE:

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.cm as cm
import plotting as myplt
from scipy.optimize import curve_fit
from tabulate import tabulate

system = sys.argv[1]
comp_cells = sys.argv[2]
buffer = sys.argv[3]

flush = sys.stdout.flush

# ----------------------------------------
# FUNCTIONS:

# print a string
def ffprint(string):
    print('%s' %(string))
    flush()

# compute residual sum of squares
def rss(y1,y2):
    return np.sum(np.power(y2-y1,2))

# compute residual sum of squares
def r2(ydata,yfit):
    return 1-rss(ydata,yfit)/np.sum(np.power(ydata-np.mean(ydata),2))

# define michaelis menten equation to get the kcat and Km with error
def mm(s,Vmax,Km):
    return Vmax*s/(Km+s)

# define michaelis menten equation to get the kcat/Km with error
def mm_kcat_Km(s,a,b):
    return a*s/(1.0+(s/b))

# define exponential with offset
def fit_exp(t,a,kobs,c):
    return a*(1.0-np.exp(-kobs*t))+c

# define exponential with offset
def fit_exp_linear(t,a,kobs,c,d):
    return a*(1.0-np.exp(-kobs*t))+ c + d*t

# define exponential with offset
def fit_bi_exp(t,a,kobs1,kobs2,c):
    return a*(1.0-np.exp(-kobs1*t)-np.exp(-kobs2*t))+ c

# ----------------------------------------
# MAIN PROGRAM:

# Loading in the datafiles
helicaseTrial01 = np.loadtxt('rep1_%s_%s_%s_data_only.dat' %(system,comp_cells,buffer))
helicaseTrial02 = np.loadtxt('rep2_%s_%s_%s_data_only.dat' %(system,comp_cells,buffer))
helicaseTrial03 = np.loadtxt('rep3_%s_%s_%s_data_only.dat' %(system,comp_cells,buffer))

# Determining the shape of the datafiles
rows = helicaseTrial01.shape[0]
columns = helicaseTrial01.shape[1]

# Creating the Trials
helicaseTrials = np.empty((3,rows,columns))
helicaseTrials[0,:,:] = helicaseTrial01
helicaseTrials[1,:,:] = helicaseTrial02
helicaseTrials[2,:,:] = helicaseTrial03

# Determining the Time for the trials
nSteps = len(helicaseTrials[0,:,0])
time = np.zeros(nSteps)
for i in range(nSteps):
    time[i] = i

# Plotting Parameters for each plot created
fig1 = plt.figure(figsize=[11,9], tight_layout={'pad':0}, facecolor='w', edgecolor='k')
fig1, ax1 = plt.subplots()
ax1.grid(b=True, which='both', axis='both', color='#808080', linestyle='--')
ax1.set_xlabel(xlabel='Time (s)',size=20)
ax1.set_ylabel(ylabel='RFU',size=20)
ax1.set_xticks(range(0,350,50))
ax1.set_yticks(range(4000,12000,1000))
ax1.tick_params(axis='both',labelsize=20)
plt.tight_layout()

fig2 = plt.figure(figsize=[11,9], tight_layout={'pad':0}, facecolor='w', edgecolor='k')
fig2, ax2 = plt.subplots()
ax2.grid(b=True, which='both', axis='both', color='#808080', linestyle='--')
ax2.set_xlabel(xlabel='Time (s)',size=20)
ax2.set_ylabel(ylabel='RFU',size=20)
ax2.set_xticks(range(0,350,50))
ax2.set_yticks(range(4000,12000,1000))
ax2.tick_params(axis='both',labelsize=20)
plt.tight_layout()

fig3 = plt.figure(figsize=[11,9], tight_layout={'pad':0}, facecolor='w', edgecolor='k')
fig3, ax3 = plt.subplots()
ax3.grid(b=True, which='both', axis='both', color='#808080', linestyle='--')
ax3.set_xlabel(xlabel='Time (s)',size=20)
ax3.set_ylabel(ylabel='RFU',size=20)
ax3.set_xticks(range(0,350,50))
ax3.set_yticks(range(4000,12000,1000))
ax3.tick_params(axis='both',labelsize=20)
plt.tight_layout()

# Defining Terms and Starting the fitting Protocol
p0_exp_linear = [5000,0.035,5000,-1]   # amplitude, kobs, offset, linear
paramNames = ["A", "kobs", "Offset", "Linear", "R2"]
trial01 = ["","","0","",""]
trial02 = ["","","1","",""]
trial03 = ["","","2","",""]
nTrials = helicaseTrials.shape[0]
nConcentrations = helicaseTrials.shape[2]
kobs = np.empty((nTrials,nConcentrations,2),dtype=float)
f1 = open('initial_params_table_%s_%s_%s_all_reps.dat' %(system,comp_cells,buffer),'w')
for trial in range(nTrials):
    for i in range(nConcentrations):
        if trial==0 and i==1:     #  wt_bl21   t407a_bl21   t407c_bl21   s411a_bl21   s411c_bl21   dm_bl21   wt_t7_oxi   wt_t7_red    dm_t7_oxi     dm_t7_red   
            start = 0             #     0          0             0         5(i=1)       5(i=1)    5(i=0,i=1)     0        62(i=1)      6(i=1)     34(i=0),24(i=1)     
            end = 300             #    300        300         202(i=1)       300          300        300        300         300          300           300    
        elif trial==1 and i==1:   #  wt_bl21   t407a_bl21   t407c_bl21   s411a_bl21   s411c_bl21   dm_bl21   wt_t7_oxi   wt_t7_red    dm_t7_oxi     dm_t7_red 
            start = 10            #  10(i=1)       0             0            0         5(i=1)        0          0         3(i=1)   5(i=0),10(i=1)    2(i=1)     
            end = 300             #    300        300           300          300          300        300        300         300          300           300    
        elif trial==2 and i==1:   #  wt_bl21   t407a_bl21   t407c_bl21   s411a_bl21   s411c_bl21   dm_bl21   wt_t7_oxi   wt_t7_red    dm_t7_oxi     dm_t7_red 
            start = 0             #     0          0          13(i=1)       5(i=1)       5(i=1)       0          0         2(i=1)         0           2(i=1)     
            end = 300             #    300        300           300          300        265(i=1)     300        300         300          300           300    
        else:
            start = 0
            end = 300
        params, var = curve_fit(fit_exp_linear,time[start:end],helicaseTrials[trial,start:end,i],p0=p0_exp_linear)
        kobs[trial,i,0] = params[1]
        kobs[trial,i,1] = var[1,1]
        R2 = r2(helicaseTrials[trial,:,i],fit_exp_linear(time[:],params[0],params[1],params[2],params[3]))
        paramList = list(params)
        paramList.append(R2)
        errorList = list(np.sqrt(np.diag(var)))
        errorList.append(0)
        if trial==0:
            table = list(zip(trial01,paramNames,paramList,errorList))
            color = next(ax1._get_lines.prop_cycler)['color']
            l1, = ax1.plot(time[:],helicaseTrials[0,:,i],'o',color=color,ms=3)
            l2, = ax1.plot(time[:],fit_exp_linear(time[:],params[0],params[1],params[2],params[3]),'--',color=color,lw=1)
            fig1.savefig('fit_exp_linear_%s_%s_%s_rep1.png' %(system,comp_cells,buffer),transparent=True,dpi=300)
        if trial==1:
            table = list(zip(trial02,paramNames,paramList,errorList))
            color = next(ax2._get_lines.prop_cycler)['color']
            l1, = ax2.plot(time[:],helicaseTrials[1,:,i],'o',color=color,ms=3)
            l2, = ax2.plot(time[:],fit_exp_linear(time[:],params[0],params[1],params[2],params[3]),'--',color=color,lw=1)
            fig2.savefig('fit_exp_linear_%s_%s_%s_rep2.png' %(system,comp_cells,buffer),transparent=True,dpi=300)
        if trial==2:
            table = list(zip(trial03,paramNames,paramList,errorList))
            color = next(ax3._get_lines.prop_cycler)['color']
            l1, = ax3.plot(time[:],helicaseTrials[2,:,i],'o',color=color,ms=3)
            l2, = ax3.plot(time[:],fit_exp_linear(time[:],params[0],params[1],params[2],params[3]),'--',color=color,lw=1)
            fig3.savefig('fit_exp_linear_%s_%s_%s_rep3.png' %(system,comp_cells,buffer),transparent=True,dpi=300)
        f1.write(tabulate(table))
    plt.close(fig1)
    plt.close(fig2)
    plt.close(fig3)
f1.close()

# Determining the kcat and KM values for the three trials
atpConc = np.array([0,1.0,10.0,30.0,100.0,200.0,400.0,600.0,800.0,1000.0])
ax = myplt.define_figure(xlabel="[ATP] ($\mu$M)",ylabel="k$_{obs}$")
ax.errorbar(atpConc[1:],kobs[0,1:,0],yerr=kobs[0,1:,1],fmt='o')
ax.errorbar(atpConc[1:],kobs[1,1:,0],yerr=kobs[1,1:,1],fmt='x')
ax.errorbar(atpConc[1:],kobs[2,1:,0],yerr=kobs[2,1:,1],fmt='v')
kcatKmNames = ["kcat", "KM","kcat/KM"]
guess = [0.05,20.0]
fitX = np.concatenate((atpConc[1:8],atpConc[1:8],atpConc[1:8]))
fitY = np.concatenate((kobs[0,1:8,0],kobs[1,1:8,0],kobs[2,1:8,0]))
params0, var0 = curve_fit(mm_kcat_Km,fitX,fitY,p0=guess)
params1, var1 = curve_fit(mm,fitX,fitY,p0=guess)
print(params0,np.sqrt(np.diag(var0)))
print(params1,np.sqrt(np.diag(var1)))
#kcatKm = params1[0]/params1[1]
#kcatKmError = np.sqrt((np.power(np.sqrt(var1[0,0])/params1[0],2))+(np.power(np.sqrt(var1[1,1])/params1[1],2)))*kcatKm
paramList = list(params1)
paramList.append(params0[0])
errorList = list(np.sqrt(np.diag(var1)))
errorList.append(np.sqrt(var0[0,0]))
table = list(zip(kcatKmNames,paramList,errorList))
f2 = open('kcat_Km_params_%s_%s_%s.dat' %(system,comp_cells,buffer),'w')
f2.write(tabulate(table))
f2.close()
x = np.arange(0,600,5)
ax.plot(x,mm_kcat_Km(x,params0[0],params0[1]))
#ax.plot(x,mm(x,params1[0],params1[1]))
plt.savefig('kobs_v_ATP_conc_for_%s_%s_%s_reps.png' %(system,comp_cells,buffer),transparent=True,dpi=300)
plt.close()
