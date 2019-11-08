import numpy as np
import MDAnalysis as md
import reach_routines as reach

#selection = "resid 1:451"
selection = "resid 1:451 or name P"

# setup MDAnalysis universe and selection
fileList = []
fileBase = "../traj/Truncated/production." 
for i in range(21,101):
    fileName = fileBase + str(i) + "/production." + str(i) + ".dcd"
    fileList.append(fileName)
coord = md.Universe("../traj/truncated.prmtop",fileList)
print("Number of frames:",coord.trajectory.n_frames)
sel = coord.select_atoms(selection)
avgPos, alignedPos  = reach.iterative_align_average_com(coord,sel,thresh=1E-8)
np.savetxt("wt_ssrna_atp_1_ssRNA_avg_structure.dat",avgPos)


#ca_sel = coord.select_atoms("name CA")
ca_sel = coord.select_atoms("name CA or name P")
covar = np.zeros((3*ca_sel.n_atoms,3*ca_sel.n_atoms),dtype=np.float64)
ca_sel.positions = avgPos
ca_sel.write("wt_ssrna_atp_1_ssRNA_avg_structure.pdb")

# loop through trajectory and compute covariance 
for ts in coord.trajectory:
    covar += np.dot(alignedPos[ts.frame-1,:,:].reshape(3*ca_sel.n_atoms,1),alignedPos[ts.frame-1,:,:].reshape(1,3*ca_sel.n_atoms))
# finish covariance
covar /= coord.trajectory.n_frames
covar -= np.dot(avgPos.reshape(3*ca_sel.n_atoms,1),avgPos.reshape(1,3*ca_sel.n_atoms))
np.savetxt("wt_ssrna_atp_1_ssRNA_com_covar.dat",covar)

