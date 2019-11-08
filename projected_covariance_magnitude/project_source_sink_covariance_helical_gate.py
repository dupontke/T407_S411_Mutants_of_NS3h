
from tabulate import tabulate
#print(tabulate([['Alice', 24], ['Bob', 19]], headers=['Name', 'Age']))
import numpy as np
from scipy import stats

# Ser364, Ile365, Lys366
sinks = [197-1,198-1,199-1]  # Ser364, Ile365, Lys366
sinks = [197-1,199-1]  # Ser364, Lys366
sinkIndeces = []
for sink in sinks:
    for i in range(3):
        sinkIndeces.append(sink*3+i)
sources = [247-1]  # Gly414
sourceIndeces = []
for source in sources:
    for i in range(3):
        sourceIndeces.append(source*3+i)
vectorIndeces = [199-1,436-1]
helixIndeces = np.arange(197,210,1).astype(int)
systems = ["wt","t407a","t407c","s411a","s411c","t407c_s411c_red","t407c_s411c_ox"]
gate1 = [197-1, 198-1, 199-1]
gate2 = [436-1, 437-1, 438-1]

covarTable = []
covarProjectionTable = []
for j, system in enumerate(systems):
    covarNorm = []
    covarProjection = []
    covarTable.append([])
    covarTable[j].append(system)
    covarProjectionTable.append([])
    covarProjectionTable[j].append(system)
    for i in range(1,4):
        # read in covariance matrix
        covarFile = system + "_ssrna_atp_" + str(i) + "_ssRNA_com_covar.dat"
        covar = np.loadtxt(covarFile)
        # read in average structure
        avgStructureFile = system + "_ssrna_atp_" + str(i) + "_ssRNA_avg_structure.dat"
        avgStructure = np.loadtxt(avgStructureFile)
        # compute p3-p1 vector
        vec = avgStructure[vectorIndeces[1],:] - avgStructure[vectorIndeces[0],:]
        vec /= np.linalg.norm(vec)   # normalize
        # gate COG vector
        vec = np.mean(avgStructure[gate2,:],axis=0) - np.mean(avgStructure[gate1,:],axis=0)
        vec /= np.linalg.norm(vec)
        # determine helical axis
        helix = avgStructure[helixIndeces,:]
        avg = np.mean(helix,axis=0)
        for l in range(helix.shape[0]):
            helix[l,:] -= avg
        e,v = np.linalg.eigh(np.dot(helix.T,helix))
        # compute norm of total covariance and append to trial list
        covarNorm.append(np.linalg.norm(covar[np.ix_(sinkIndeces,sourceIndeces)]))
        # project covariance along vector and compute norm
        #projection = []
        for source in sources:
            for sink in sinks:
                #projection.append(np.dot(covar[source*3:source*3+3,sink*3:sink*3+3],v[:,2]))
                covarProjection.append(np.linalg.norm(np.dot(covar[source*3:source*3+3,sink*3:sink*3+3],vec)))
        #covarProjection.append(np.linalg.norm(projection))
    #covarNorm = np.asarray(covarNorm)
    covarTable[j].append(np.mean(covarNorm))
    covarTable[j].append(stats.sem(covarNorm))
    covarTable[j].append(np.mean(covarNorm)/covarTable[0][1])
    covarProjectionTable[j].append(np.mean(covarProjection))
    covarProjectionTable[j].append(stats.sem(covarProjection))
    covarProjectionTable[j].append(np.mean(covarProjection)/covarProjectionTable[0][1])
print("TOTAL COVARIANCE BETWEEN SOURCES AND SINKS:")
print(tabulate(covarTable, headers=['System', 'Mean', 'Stderr', 'fold difference']))
print("\n\nPROJECTED COVARIANCE BETWEEN SOURCES AND SINKS:")
print(tabulate(covarProjectionTable, headers=['System', 'Mean', 'Stderr', 'fold difference']))
