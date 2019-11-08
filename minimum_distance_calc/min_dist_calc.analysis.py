#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ./com.analysis.py pdb_file trajectory_location start end system_descriptor

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import math
from numpy.linalg import *
import MDAnalysis
import MDAnalysis.analysis.distances
from distance_functions import *
from sel_list import *

# ----------------------------------------
# VARIABLE DECLARATION:
 
pdb = sys.argv[1]                   # point to a pdb or prmtop or psf file (untested for both prmtop and psf files)
traj_loc = sys.argv[2]              # point to the location of the trajectory files
start = int(sys.argv[3])            # integer describing
end = int(sys.argv[4])
system = sys.argv[5]

zeros = np.zeros
square = np.square
sqrt = np.sqrt
flush = sys.stdout.flush

important = 'protein or nucleic or resname A5 or resname A3 or resname U5 or resname ATP or resname MG'
nSel = len(sel)

# ----------------------------------------
# FUNCTIONS:

def ffprint(string):
    print '%s' %(string)
    flush()

# ----------------------------------------
# MAIN PROGRAM:

u = MDAnalysis.Universe(pdb)
u_important = u.select_atoms(important)

nRes = len(u_important.residues)
ffprint(nRes)

# make an atom selection to compute COM distance
u_sel = []
for i in range(nSel):
    selection0 = sel[i][0]
    selection1 = sel[i][1]
    u_sel.append([u.select_atoms(selection0),u.select_atoms(selection1)])
    ffprint('%s atoms found in %s' %(u_sel[i][0].n_atoms, u_sel[i][0].residues))
    ffprint('%s atoms found in %s' %(u_sel[i][1].n_atoms, u_sel[i][1].residues))

hist_min = 0
hist_max = 100
bin_size = 0.1
num_bins = int((hist_max-hist_min)/bin_size)
dist_hist = np.zeros(num_bins,dtype=float)
box = u.dimensions[:3]
ffprint('the box dimensions are %s' %(box))

out1 = open('%03d.%03d.%s.min_dist_calc.dat' %(int(sys.argv[3]),end,system),'w')

count = 0
nSteps = 0
while start <= end:
    ffprint('Loading trajectory %s' %(start))
    u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
    nSteps += len(u.trajectory)

    for ts in u.trajectory:
        if ts.frame%1000 == 0:
            ffprint('Working on timestep %d of trajectory %d' %(ts.frame, start))

        for i in range(nSel):
            res0 = u_sel[i][0].positions
            res1 = u_sel[i][1].positions
            min_dist2 = 9999.
            for j in range(len(u_sel[i][0])):
                atom0 = res0[j]
                for k in range(len(u_sel[i][1])):
                    atom1 = res1[k]
                    dist2 = euclid_dist(atom0,atom1)[0]
                    if dist2 < min_dist2:
                        min_dist2 = dist2
            #dist,dist2 = euclid_dist(com0,com1)
            #dist = math.sqrt(computePbcDist2(com0,com1,box))
            out1.write('%10.6f    ' %(min_dist2))
        out1.write('\n')
    start +=1

out1.close()
ffprint(nSteps)


# print histogram data into separate file
#with open('%03d.%03d.prob_density_hist.dat' %(int(sys.argv[3]),end),'r') as file:
#    dat = np.loadtxt(file,dtype=np.float)

#dist_bin = int((dist-hist_min)/bin_size)
#if hist_min<dist<hist_max:
#    dist_hist[dist_bin] += 1
#    for j in range(num_bins):
#        dist_hist[i] /= 4*math.pi*((i+0.5)*bin_size+hist_min)**2
#out2 = open('%03d.%03d.prob_density_hist.dat' %(int(sys.argv[3]),end),'w')
#    out2.write('%10.6f   %10.6f ' %(i*bin_size+hist_min,dist_hist[i]))
#    out2.write('\n')
#out2.close()




