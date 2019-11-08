#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:
#./nucleophilic_water_plotting.py config_file

# ----------------------------------------
# PREAMBLE:

import numpy as np
import numpy.ma as ma
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

sqrt = np.sqrt
zeros = np.zeros
nullfmt = NullFormatter()
prob_density_cmap = plt.cm.get_cmap('Blues')
prob_density_cmap.set_under('w',alpha=0)
free_energy_cmap = plt.cm.get_cmap('jet')

config_file = sys.argv[1]

k = 0.001987 # Kcal K^-1 mol^-1
T = 300. # K
kT = k*T
four_pi =4.0* np.pi

# ----------------------------------------
# FUNCTIONS

necessary_parameters =  ['system_descriptor','nucleophile_water_datafile','average_water_datafile','attack_distance_cutoff','attack_angle_cutoff','dipole_angle_cutoff','start','stop','step']
all_parameters =        ['system_descriptor','nucleophile_water_datafile','average_water_datafile','attack_distance_cutoff','attack_angle_cutoff','dipole_angle_cutoff','start','stop','step','attack_distance_bin_width','angle_bin_width','delta_write','homemade_time_range_bool','homemade_time_range_list','plotting_bool','nucleophile_water_postanalysis_output','average_water_postanalysis_output','time_chunk_output_filename']   #'Free_energy_bool',

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['delta_write'] = float(0.002)		# Units: ns
	
        parameters['plotting_bool'] = False
	#parameters['Free_energy_bool'] = False
	
        parameters['attack_distance_bin_width'] = float(0.05)		# Units: Angstroms
	parameters['angle_bin_width'] = float(1.50)		# Units: Degrees
	parameters['homemade_time_range_bool'] = False
	parameters['homemade_time_range_list'] = None
	parameters['nucleophile_water_postanalysis_output'] = 'nucleophilic_waters.dat'
	parameters['average_water_postanalysis_output'] = 'average_waters.dat'
	parameters['plot_x_range_tuple'] = (2.50,8.00)
	parameters['plot_y_range_tuple'] = (25.0,180.0)
	parameters['time_chunk_output_filename'] = 'data_table.dat'

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)

def Nucleophilic_waters(dist_cutoff,ang_cutoff,dipole_cutoff,data_array):
        """ function to determine which waters (rows) in data_array are lytic; returns the total number of nucleophilic waters found in the pocket, total number of frames that possess a nucleophilic water, and an array with the corresponding indices of waters (rows) that passed the conditions to be defined lytic
        Arguments: 
            dist_cutoff = float; the geometric cutoff used for the P_gamma-O_water distance
            angle_cutoff = float; the geometric cutoff used for the O_water-P_gamma-O_beta,gamma angle
            dipole_cutoff = float; the geometric cutoff used for the Water dipole vector- P_gamma-O_beta,gamma vector angle
            data_aray = np.array; contains all of the necessary data for this analysis; same organization as the data file read into this script 
        """
        ### APPLYING OLD CUTOFFS; ONLY GEOMETRIC CUTOFFS
        geom_nucl_water_indices = np.array([i for i in range(len(data_array)) if (data_array[i,2]<dist_cutoff and data_array[i,3] > ang_cutoff)])
        geom_wat_count = len(geom_nucl_water_indices)  # the number of elements in the zeroth element in the tuple corresponds to the total number of waters that passed the conditions
        geom_frame_count = len(set(data_array[geom_nucl_water_indices][:,0]))      # grabbing lines from data_array with the corresponding indices stored in nucl_water_indices, grabbing the zeroth column in each row of this abbreviated data_array, making a set out of this array. This set now contains the frame numbers that have waters that passed the conditions; the length of the set corresponds to the frame_count
        
        #geom_nucl_water_indices = np.where(np.logical_and(data_array[:,2]<dist_cutoff,data_array[:,3]>ang_cutoff))   # searches through the data_array for rows that pass both conditions; creates a weird tuple where the first element in the tuple corresponds to the indices of data_array where waters have passed conditions
        #geom_wat_count = len(geom_nucl_water_indices[0])  # the number of elements in the zeroth element in the tuple corresponds to the total number of waters that passed the conditions
        #geom_frame_count = len(set(data_array[geom_nucl_water_indices[0]][:,0]))      # trickly line; grabbing lines from data_array with the corresponding indices stored in nucl_water_indices, grabbing the zeroth column in each row of this abbreviated data_array, making a set out of this array. This set now contains the frame numbers that have waters that passed the conditions; the length of the set corresponds to the frame_count

        ### APPLYING NEW CUTOFFS 
        nucl_water_indices = np.array([i for i in range(len(data_array)) if (data_array[i,2] < dist_cutoff and data_array[i,3] > ang_cutoff and data_array[i,4] > dipole_cutoff)])
        wat_count = len(nucl_water_indices)
        frame_count = len(set(data_array[nucl_water_indices][:,0]))
	#temp_data = np.array([i for i in data_array[geom_nucl_water_indices[0]] if i[4] > dipole_cutoff])
        
        return geom_wat_count, geom_frame_count, geom_nucl_water_indices, wat_count, frame_count, nucl_water_indices	# RETURNING THREE OBJECTS PER CUTOFF STEPS; THE TOTAL NUMBER OF NUCLEOPHILIC WATERS FOUND IN THE POCKET; THE TOTAL NUMBER OF FRAMES WITH A NUCLEOPHILIC WATER; THE ARRAY INDICES CORRESPONDING TO ROWS IN THE DATA_ARRAY THAT PASS THE CONDITIONS

# ----------------------------------------
# MAIN
# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

start = int(parameters['start'])
stop = int(parameters['stop'])
step = int(parameters['step'])
delta_write = float(parameters['delta_write'])

attack_distance_cutoff = float(parameters['attack_distance_cutoff'])
attack_distance_minimum = float(parameters['attack_distance_minimum'])
attack_distance_maximum = float(parameters['attack_distance_maximum'])
attack_distance_bin_width = float(parameters['attack_distance_bin_width'])

attack_angle_cutoff = float(parameters['attack_angle_cutoff'])
dipole_angle_cutoff = float(parameters['dipole_angle_cutoff'])

angle_minimum = float(parameters['angle_minimum'])
angle_maximum = float(parameters['angle_maximum'])
angle_bin_width = float(parameters['angle_bin_width'])

# ----------------------------------------
# PREP TIME RANGES TO ANALYZE THE NUCLEOPHILIC/AVERAGE WATER DATA
frame_ranges = []
for i in range(start,stop,step):	# i has units of ns
	j = i + step 			# j has units of ns
	first_frame = int(i/delta_write)
	last_frame = int(j/delta_write)-1
	nSteps = len(range(first_frame,last_frame))+1
	frame_ranges.append([nSteps,i,j,first_frame,last_frame])

if parameters['homemade_time_range_bool']:
	for i in range(len(parameters['homemade_time_range_list'])): # assumes the time_range_list is a list of lists....
		first_frame = int(parameters['homemade_time_range_list'][i][0]/delta_write)
		last_frame = int(parameters['homemade_time_range_list'][i][1]/delta_write)-1
		nSteps = len(range(first_frame,last_frame))+1
		frame_ranges.append([nSteps,parameters['homemade_time_range_list'][i][0],parameters['homemade_time_range_list'][i][1],first_frame,last_frame])

# ----------------------------------------
# DETERMINE THE NUMBER OF BINS TO BE USED IN PLOTTING THE NUCLEOPHILIC DATA
num_bins = [int((attack_distance_maximum-attack_distance_minimum)/attack_distance_bin_width),int((angle_maximum-angle_minimum)/angle_bin_width)]
delta_attack_dist = (attack_distance_maximum - attack_distance_minimum)/num_bins[0]
delta_angle = (angle_maximum - angle_minimum)/num_bins[1]
print num_bins,delta_attack_dist,delta_angle

dist_edges = zeros(num_bins[0]+1,dtype=np.float64)
angle_edges = zeros(num_bins[1]+1,dtype=np.float64)
for i in range(num_bins[0]+1):
	dist_edges[i] = attack_distance_minimum + i*delta_attack_dist
for i in range(num_bins[1]+1):
	angle_edges[i] = angle_minimum + i*delta_angle

# ----------------------------------------
# LOAD THE NUCLEOPHILIC WATER DATAFILE INTO MEMORY
nucleophilic_water_data = np.loadtxt(parameters['nucleophile_water_datafile'])

# ----------------------------------------
# One dimensional histograms of all attack distance, attack angle, and dipole angle metrics
events,edges,patches = plt.hist(nucleophilic_water_data[:,2],bins=num_bins[0],histtype='bar',normed=True,edgecolor='black')
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel(r'O$_{wat}$ - P$_{\gamma}$ Distance ($\AA$)',size=20)
plt.ylabel('Probability Density',size=20)
plt.tight_layout()
plt.savefig('%s.full_data.attack_distance.png'%(parameters['system_descriptor']),dpi=600,transparent=True)
plt.close()

events,edges,patches = plt.hist(nucleophilic_water_data[:,3],bins=num_bins[1],histtype='bar',normed=True,edgecolor='black')
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel(r'O$_{wat}$-P$_{\gamma}$-O$_{\beta,\gamma}$ Angle (Deg)',size=20)
plt.ylabel('Probability Density',size=20)
plt.tight_layout()
plt.savefig('%s.full_data.attack_angle.png'%(parameters['system_descriptor']),dpi=600,transparent=True)
plt.close()

events,edges,patches = plt.hist(nucleophilic_water_data[:,4],bins=num_bins[1],histtype='bar',normed=True,edgecolor='black')
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel(r'Dipole Angle (Deg)',size=20)
plt.ylabel('Probability Density',size=20)
plt.tight_layout()
plt.savefig('%s.full_data.dipole_angle.png'%(parameters['system_descriptor']),dpi=600,transparent=True)
plt.close()

# ----------------------------------------
# BEGIN DATA ANALYSIS OF THE NUCLEOPHILIC WATER DATA
with open(parameters['nucleophile_water_postanalysis_output'],'w') as W, open(parameters['time_chunk_output_filename'],'w') as Y:       # open the necessary files to contain all results

        W.write('Distance Cutoff: Less than %f \nAttack Angle Cutoff: Greater than %f\nDipole Angle Cutoff: Greater than %f\n\nTotal Dataset\n	Smallest P_gamma-O_water distance (nucleophilic): %f\n	Largest P_gamma-O_water distance: %f\n	Smallest O_water-P_gamma-O_beta,gamma angle: %f\n	Largest O_water-P_gamma-O_beta,gamma angle (nucleophilic): %f\n Smallest dipole angle: %f\n Largest dipole angle (nucleophilic): %f\n\nNumber of bins in the distance dimension: %d\nNumber of bins in the angle dimensions: %d\n Distance bin dimension: %f; Angle bin dimensions: %f\n\n' %(attack_distance_cutoff,attack_angle_cutoff,dipole_angle_cutoff,np.min(nucleophilic_water_data[:,2]),np.max(nucleophilic_water_data[:,2]),np.min(nucleophilic_water_data[:,3]),np.max(nucleophilic_water_data[:,3]),np.min(nucleophilic_water_data[:,4]),np.max(nucleophilic_water_data[:,4]),num_bins[0],num_bins[1],delta_attack_dist,delta_angle))   # write to output file the parameters used to develop the analysis 

	for i in frame_ranges:  # loop over all of the ranges of trajectory (previously determined)
		first_index = np.argwhere(nucleophilic_water_data[:,0]==i[3])[0][0]	# SINCE THE DATA IS NOT ORGANIZED BY TIMESTEP, I NEED TO DETERMINE THE CORRECT INDEX FOR THE FIRST FRAME OF THE TIME RANGE
		last_index = np.argwhere(nucleophilic_water_data[:,0]==i[4])[-1][0]	# SIMILAR TO ABOVE; DETERMINE THE INDEX OF THE LAST FRAME IN THE TIME RANGE
		nValues = len(range(first_index,last_index))+1      # number of waters found in the active site
		W.write('Analyzing %d to %d ns; corresponds to frames %d to %d, which have indices/line numbers %d to %d\nCorresponds to %d steps\nCorresponds to %d waters\n'%(i[1],i[2],i[3],i[4],first_index,last_index,i[0],nValues))
                time_range_data = nucleophilic_water_data[first_index:last_index]

		# ----------------------------------------
		# CALCULATE THE PROBABILITY OF OBSERVING A NUCLEOPHILIC WATER IN THE TIME RANGE; OUTPUT TO FILE
		geometric_water_count, geometric_frame_count, geometric_nucl_water_indices, all_cutoffs_wat_count, all_cutoffs_frame_count, all_cutoffs_nucl_water_indices = Nucleophilic_waters(attack_distance_cutoff,attack_angle_cutoff,dipole_angle_cutoff,time_range_data)   # calculate the number of nucleophilic waters using the geometric definition and cutoffs
		
                wat_prob = all_cutoffs_wat_count/float(nValues)     # turn number into a prob.
		wat_prob_error = sqrt(all_cutoffs_wat_count)/float(nValues) # calculate the std. error of the mean assuming counting experiment error
		frame_prob = all_cutoffs_frame_count/float(i[0])     # turn number into a prob.
		frame_prob_error = sqrt(all_cutoffs_frame_count)/float(i[0])  # calculate the std. error of the mean assuming counting experiment error

		Y.write('%10d   %10d   %10d   %10f   %10f   %10d   %10f   %f\n' %(i[1],i[2],all_cutoffs_wat_count,wat_prob,wat_prob_error,all_cutoffs_frame_count,frame_prob,frame_prob_error))
		W.write('Total number of nucleophilic waters in this time range: %d\nProbability of waters being nucleophilic (includes frames with multiple nucleophilic waters): %f\nProbability error arrived at by sqrt rule of counting experiments: %f\nTotal number of frames with nucleophilic waters: %d\nProbability of a frame having a nucleophilc water: %f\nProbability error arrived at by sqrt rule of counting experiments: %f\n\n' %(all_cutoffs_wat_count,wat_prob,wat_prob_error,all_cutoffs_frame_count,frame_prob,frame_prob_error))  # print out to global file the results for this time range.

		# ----------------------------------------
		# PLOT THE PROB DENSITIES OF THE NUCLEOPHILIC WATER DATA; 
		if parameters['plotting_bool']:
                        os.mkdir('%04d.%04d.figures'%(i[1],i[2]))
                        # PLOT attack distance
                        events,edges,patches = plt.hist(time_range_data[all_cutoffs_nucl_water_indices][:,2],bins=num_bins[0],range=(attack_distance_minimum,attack_distance_maximum),normed=True,edgecolor='black')
                        plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
                        plt.title('Time Range: %s to %s'%(i[1],i[2]),size=20)
                        plt.xlabel(r'O$_{wat}$ - P$_{\gamma}$ Distance ($\AA$)',size=20)
                        plt.ylabel('Probability Density',size=20)
                        plt.xlim((0,10))
                        plt.ylim((0,3.0))
			plt.tight_layout()
                        plt.savefig('%04d.%04d.figures/%04d.%04d.attack_distance.png'%(i[1],i[2],i[1],i[2]),dpi=600,transparent=True)
                        plt.close()

                        # PLOT attack angle
                        events,edges,patches = plt.hist(time_range_data[all_cutoffs_nucl_water_indices][:,3],bins=num_bins[1],range=(angle_minimum,angle_maximum),normed=True,edgecolor='black')
                        plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
                        plt.title('Time Range: %s to %s'%(i[1],i[2]),size=20)
                        plt.xlabel(r'O$_{wat}$-P$_{\gamma}$-O$_{\beta,\gamma}$ Angle (Deg)',size=20)
                        plt.ylabel('Probability Density',size=20)
                        plt.xlim((0,180))
                        plt.ylim((0,0.15))
			plt.tight_layout()
                        plt.savefig('%04d.%04d.figures/%04d.%04d.attack_angle.png'%(i[1],i[2],i[1],i[2]),dpi=600,transparent=True)
                        plt.close()

                        # PLOT dipole angle
                        events,edges,patches = plt.hist(time_range_data[all_cutoffs_nucl_water_indices][:,4],bins=num_bins[1],range=(angle_minimum,angle_maximum),normed=True,edgecolor='black')
                        plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
                        plt.title('Time Range: %s to %s'%(i[1],i[2]),size=20)
                        plt.xlabel(r'Dipole Angle (Deg)',size=20)
                        plt.ylabel('Probability Density',size=20)
                        plt.xlim((0,180))
                        plt.ylim((0,0.05))
			plt.tight_layout()
                        plt.savefig('%04d.%04d.figures/%04d.%04d.dipole_angle.png'%(i[1],i[2],i[1],i[2]),dpi=600,transparent=True)
                        plt.close()

                        # PLOT 2d histogram of collective variables; NOTE: HIST2D FUNCTION CALCULATES THE PROB DENSITY AS NORMALIZED BY THE TOTAL NUMBER OF DATA POINTS (INSTEAD OF NSTEPS) * DX * DY
                        # atk angle vs atk dist
			counts,xedges,yedges,image = plt.hist2d(nucleophilic_water_data[first_index:last_index,2],nucleophilic_water_data[first_index:last_index,3],bins=[num_bins[0],num_bins[1]],range=[[attack_distance_minimum,attack_distance_maximum],[angle_minimum,angle_maximum]],cmap=prob_density_cmap,normed=True,vmin=0.001,vmax=0.020)	# 2d histogram (aka heat map) of the geometric collective variables; probability density 
			cb1 = plt.colorbar(extend='max')		# extend='max'
			cb1.set_label('Probability Density',size=20)
			plt.xlabel(r'P$_{\gamma}$ - O$_{wat}$ Distance ($\AA$)',size=20)
			plt.ylabel(r'O$_{\beta,\gamma}$-P$_{\gamma}$-O$_{wat}$ Angle (Deg)',size=20)
			plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
			plt.xlim(parameters['plot_x_range_tuple'])
			plt.ylim(parameters['plot_y_range_tuple'])
			plt.tight_layout()
			plt.savefig('%04d.%04d.figures/%04d.%04d.atk_angle_atk_dist.hist2d.png'%(i[1],i[2],i[1],i[2]),transparent=True,dpi=600)
			plt.close()

                        # dipole angle vs atk dist
			counts,xedges,yedges,image = plt.hist2d(nucleophilic_water_data[first_index:last_index,2],nucleophilic_water_data[first_index:last_index,4],bins=[num_bins[0],num_bins[1]],range=[[attack_distance_minimum,attack_distance_maximum],[angle_minimum,angle_maximum]],cmap=prob_density_cmap,normed=True,vmin=0.001,vmax=0.020)	# 2d histogram (aka heat map) of the geometric collective variables; probability density 
			cb1 = plt.colorbar(extend='max')		# extend='max'
			cb1.set_label('Probability Density',size=20)
			plt.xlabel(r'P$_{\gamma}$ - O$_{wat}$ Distance ($\AA$)',size=20)
			plt.ylabel(r'Dipole Angle (Deg)',size=20)
			plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
			plt.xlim(parameters['plot_x_range_tuple'])
			plt.ylim(parameters['plot_y_range_tuple'])
			plt.tight_layout()
			plt.savefig('%04d.%04d.figures/%04d.%04d.dipole_angle_atk_dist.hist2d.png'%(i[1],i[2],i[1],i[2]),transparent=True,dpi=600)
			plt.close()

                        # dipole angle vs atk dist
			counts,xedges,yedges,image = plt.hist2d(nucleophilic_water_data[first_index:last_index,3],nucleophilic_water_data[first_index:last_index,4],bins=[num_bins[1],num_bins[1]],range=[[angle_minimum,angle_maximum],[angle_minimum,angle_maximum]],cmap=prob_density_cmap,normed=True,vmin=0.00001,vmax=0.00020)	# 2d histogram (aka heat map) of the geometric collective variables; probability density 
			cb1 = plt.colorbar(extend='max')		# extend='max'
			cb1.set_label('Probability Density',size=20)
			plt.xlabel(r'O$_{\beta,\gamma}$-P$_{\gamma}$-O$_{wat}$ Angle (Deg)',size=20)
			plt.ylabel(r'Dipole Angle (Deg)',size=20)
			plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
			plt.xlim(parameters['plot_y_range_tuple'])
			plt.ylim(parameters['plot_y_range_tuple'])
			plt.tight_layout()
			plt.savefig('%04d.%04d.figures/%04d.%04d.dipole_angle_atk_angle.hist2d.png'%(i[1],i[2],i[1],i[2]),transparent=True,dpi=600)
			plt.close()

		## ----------------------------------------
		## CALCULATE, WRITE TO OUTPUT, AND PLOT THE RELATIVE FREE ENERGY OF THE NUCLEOPHILIC WATER DATA; NORMALIZED FOR SPHERICAL VOLUME
		#if parameters['Free_energy_bool']:
		#	# DEFINE CONSTANTS
		#	fe_divisor = nValues*delta_attack_dist*delta_angle*four_pi
		#	# ARRAY DECLARATION
		#	counts = np.zeros(tuple(num_bins),dtype=np.float64)
		#	fe_counts = np.zeros(tuple(num_bins),dtype=np.float64)
		#	# LOOP THROUGH AND HISTOGRAM DATA
		#	for a in nucleophilic_water_data[first_index:last_index]:
		#		attack_distance_index = int((a[2] - attack_distance_minimum)/delta_attack_dist)
		#		attack_angle_index = int((a[3] - attack_angle_minimum)/delta_attack_angle)
		#		
		#		if attack_distance_index not in range(num_bins[0]+1) or attack_angle_index not in range(num_bins[1]+1):
		#			print 'A specific water is trying to be binned in a bad bin. Here is the data that is bad: ', i, attack_distance_index
		#			sys.exit()
		#		if attack_distance_index == num_bins[0]: 
		#			attack_distance_index = -1 
		#		if attack_angle_index == num_bins[1]: 
		#			attack_angle_index = -1 
		#		else:
		#			counts[attack_distance_index][attack_angle_index] += 1
		#			fe_counts[attack_distance_index][attack_angle_index] += 1/(fe_divisor*a[2]**2)	# nSteps*delta_attack_dist*delta_attack_angle normalizes the probability density; 4.0*pi*a[2]**2 volume corrects the spherical volume these waters were found in. 
		#	# FROM THE HISTOGRAMMED ARRAYS, CALCULATE THE FREE ENERGY SURFACE
		#	for j in range(num_bins[0]):
		#		for k in range(num_bins[1]):
		#				fe_counts[j][k] = -kT*np.log(fe_counts[j][k])
		#	fe_counts -= np.ndarray.min(fe_counts)		# Make lowest energy bin the zero point reference
		#	
		#	# PCOLORMESH ERRORS WHEN THE C MATRIX HAS NANS OR INFS; CREATE A MASKED ARRAY SO THAT PCOLORMESH WILL IGNORE INFS IN THE FE_COUNTS ARRAY
		#	masked_fe_counts = ma.masked_where(np.isinf(fe_counts),fe_counts)
		#	# PLOT FREE ENERGY SURFACE
		#	free_plot = plt.pcolormesh(dist_edges,angle_edges,masked_fe_counts.T,cmap=free_energy_cmap)
		#	cb1 = plt.colorbar(extend='max')		# extend='max'
		#	cb1.set_label(r'Relative Free Energy (kcal mol$^{-1}$)')
		#	plt.xlabel(r'P$_{\gamma}$ - O$_{water}$ Distance ($\AA$)',size=20)
		#	plt.ylabel(r'O$_{\beta,\gamma}$-P$_{\gamma}$-O$_{water}$ Angle (Deg)',size=20)
		#	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
		#	plt.xlim(parameters['plot_x_range_tuple'])
		#	plt.ylim(parameters['plot_y_range_tuple'])
		#	plt.tight_layout()
		#	plt.savefig('%04d.%04d.nucleophilic_wat.free_energy.png'%(i[1],i[2]),transparent=True,dpi=600)
		#	plt.close()

		#	with open('%04d.%04d.nucleophilic_wat.free_energy.dat'%(i[1],i[2]),'w') as Y:
		#		np.savetxt(Y,fe_counts)
		
# ----------------------------------------
# CALCULATE THE AVERAGE (AND STD DEV) NUMBER OF WATERS WITHIN THE BINDING POCKET
avg_water_data = np.loadtxt(parameters['average_water_datafile'])
with open(parameters['average_water_postanalysis_output'],'w') as W:
	for i in frame_ranges:
		first_index = i[3]
		last_index = i[4]
		avg_waters = np.sum(avg_water_data[first_index:last_index])
		avg_waters /= i[0]
		std_waters = np.std(avg_water_data[first_index:last_index])
		W.write('Time range: %d to %d ns; average number of waters in the ATP binding pocket: avg %f st dev %f \n' %(i[1],i[2],avg_waters,std_waters))
	
