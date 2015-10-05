# -*- coding: utf-8 -*-
"""
Created on Wed Jul 08 19:42:24 2015

@author: dan
"""

import numpy as np
import gzip
import sys

def addToCounts(vec, counts, minAC, maxAC):
    indices = [i for i,x in enumerate(vec) if x=='1']
    mac = len(indices)
    if mac<minAC or mac>maxAC:
        return
    global lineCount
    lineCount += 1
    for i in indices:
        counts[i] += 1
    
    
lineCount = 0
if __name__ == "__main__":
    chrNumber = sys.argv[1]
    minAC = int(sys.argv[2])
    maxAC = int(sys.argv[3])
    print chrNumber
    counts = [0]*5008
    filename = './data/1000GP_Phase3_chr'+str(chrNumber)+'.hap.gz'
    print filename
    index = 0
    with gzip.open(filename) as f:
        for line in f:
            line = line.replace(" ","").rstrip()
            addToCounts(line, counts, minAC, maxAC)
    print "FINISHED"
    np.savetxt('./counts_' + str(minAC) + '_' + str(maxAC) + '/counts_chr'+str(chrNumber)+'.csv', counts, delimiter=",")
