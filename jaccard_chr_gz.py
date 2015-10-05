# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 19:42:24 2015

@author: dan
"""

import numpy as np
import gzip
import sys

def addToJaccardMatrix(vec, U, V, counts, minAC, maxAC):
    indices = [i for i,x in enumerate(vec) if x=='1']
    mac = len(indices)
    if mac<minAC or mac>maxAC:
        return
    for i in indices:
        counts[i] += 1
    for i in indices:
        for j in indices:
            V[i,j] += 1
        U[:,i] += 1
        U[i,:] += 1
    
if __name__ == "__main__":
    chrNumber = 22
    minAC = int(10)
    maxAC = int(10)
    print chrNumber
    counts = [0]*5008
    U = np.zeros((5008,5008))
    V = np.zeros((5008,5008))
    filename = './data/1000GP_Phase3_chr'+str(chrNumber)+'.hap.gz'
    print filename
    with gzip.open(filename) as f:
        for line in f:
            line = line.replace(" ","").rstrip()
            addToJaccardMatrix(line, U, V, counts, minAC, maxAC)
    U = U-V
    np.savetxt('./output_' + str(minAC) + '_' + str(maxAC) + '/chr'+str(chrNumber)+'_U.csv.gz', U, delimiter=",")
    np.savetxt('./output_' + str(minAC) + '_' + str(maxAC) + '/chr'+str(chrNumber)+'_V.csv.gz', V, delimiter=",")
    np.savetxt('./output_' + str(minAC) + '_' + str(maxAC) + '/counts_chr'+str(chrNumber)+'.csv', counts, delimiter=",")

