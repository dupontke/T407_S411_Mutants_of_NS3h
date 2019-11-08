#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from plotting_functions import *
from sel_list import *

zeros = np.zeros

dat = sys.argv[1]  
system = sys.argv[2]

nSel = len(sel)

bin_size = 0.01
k = 0.001987 # Kcal K^-1 mol^-1
T = 300. # K
kT = k*T
boltz = 2*kT
four_pi = 4*np.pi


# ----------------------------------------
# MAIN PROGRAM:
# ----------------------------------------

# Load in data_file into a numpy array
datalist = np.loadtxt(dat)

nSteps = len(datalist[:,0])          
print 'Number of selections: %d, number of steps: %d' %(nSel,nSteps)

time = np.zeros(nSteps)
for i in range(nSteps):
	time[i] = i*0.002		# units of time in ns; each frame is separated by 0.002 ns 

rows = datalist.shape[0]
columns = datalist.shape[1]
if columns != nSel:
        sys.exit()
print "Number of data columns:", columns
print "Number of data rows:", rows


for i in range(nSel):
	selection = sel[i][2]
	scat_hist(time[:],datalist[:,i],'k','Time (ns)','Minimum Distance','%02d.%s.%s' %(i,selection,system),'minimum_distance',yunits='$\AA$')
	hist1d(datalist[:,i],'Minimum Distance','%02d.%s.%s' %(i,selection,system),'minimum_distance',norm=True,xunits='$\AA$')
	
	# Loop through each column and create a histogram and the probability density
	out1 = open('%02d.%s.%s.minimum_distance_calc.prob_density_hist.dat' %(i,selection,system),'w')
	
	# determine domain of data
	max_val = np.amax(datalist[:,i])
	min_val = np.amin(datalist[:,i])
	num_bins = int((max_val-min_val)/bin_size)+1

	# allocate probability and x2 arrays
	prob = zeros((num_bins),dtype=np.float)
	half_bin = zeros((num_bins),dtype=np.float)

        # create histogram of data
	for j in range(rows):
		current_bin = int((datalist[j,i] - min_val)/bin_size)
		prob[current_bin] += 1.

	# finish probability density
	prob /= bin_size*rows

	# obtain x_axis values for scatter plot
	for j in range(num_bins):
		half_bin[j] = min_val+(j*bin_size)+(bin_size/2)

	# write out to prob_density_hist.dat file
	for j in range(num_bins):
                out1.write('%10.6f \n' %(prob[j]))
        out1.close()

        # scatter plot for probability density
	plot_1d(half_bin,prob,'k','Minimum Distance','Probability Density','%02d.%s.%s' %(i,selection,system),'minimum_distance_prob_density',xunits='$\AA$')
	hist2d(datalist[:,1],datalist[:,2],'R50 gtp_gamma_phosphate Minimum Distance','S143 gtp_alpha_phosphate Minimum Distance',500,'R50gtp_gamma_S143gtp_alpha.%s' %(system),'joint_probability_density',norm=True,xunits='$\AA$',yunits='$\AA$')
	hist2d(datalist[:,1],datalist[:,3],'R50 gtp_gamma_phosphate Minimum Distance','S143 gtp_beta_phosphate Minimum Distance',500,'R50gtp_gamma_S143gtp_beta.%s' %(system),'joint_probability_density',norm=True,xunits='$\AA$',yunits='$\AA$')
