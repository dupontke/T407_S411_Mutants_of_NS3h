
from tabulate import tabulate
import numpy as np
from scipy import stats

# Lys388 is 221-1
sinks = [58-1,220-1] # arginines Arg225, Arg387
sinks = [220-1]
sinkIndeces = []
for sink in sinks:
    for i in range(3):
        sinkIndeces.append(sink*3+i)
sources = [118-1,247-1, 289-1, 293-1]  # Glu285, Gly414, Gln456, Arg460
sources = [247-1]
sourceIndeces = []
for source in sources:
    for i in range(3):
        sourceIndeces.append(source*3+i)
phosphates = [452,454]

systems = ["wt","t407a","t407c","s411a","s411c","t407c_s411c_red","t407c_s411c_ox"]

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
        p3_p1 = avgStructure[phosphates[1],:] - avgStructure[phosphates[0],:]
        p3_p1 /= np.linalg.norm(p3_p1)   # normalize
        # compute norm of total covariance and append to trial list
        covarNorm.append(np.linalg.norm(covar[np.ix_(sinkIndeces,sourceIndeces)]))
        # project covariance along p3_p1 vector and compute norm
        for source in sources:
            for sink in sinks:
                covarProjection.append(np.linalg.norm(np.dot(covar[source*3:source*3+3,sink*3:sink*3+3],p3_p1)))
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
