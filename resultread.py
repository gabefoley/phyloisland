# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 09:39:05 2018

@author: Rowan
"""

import os
import glob
from Bio import SearchIO
import csv

paths = ['DDQSY', 'FJWHO', 'ITMTT']

def getStartPostion(record, hit_start):
    print ('Lets get this start position')
    print (record)
    print (record.sequence[hit_start])
    print (record.sequence[hit_start: hit_start + 100])

def resultRead(paths, record=None):
    hmm_dict = {}
    for path in paths:
        if glob.glob(os.path.join(path, '*.csv')):
            break
        for infile in glob.glob( os.path.join(path, '*.fasta')):
            qresult = SearchIO.read(infile, 'hmmer3-text')
            for i in range(len(qresult.hsps)):
                    try:
                        if "forward" in infile:
                            hsp = qresult[0][i]
                            hmm_dict[infile[6:-24] + "forward" + '_' +str(i)] = str(3*hsp.env_start) +':'+ str(3*hsp.env_end)
                        else:
                            hsp = qresult[0][i]
                            hmm_dict[infile[6:-24] + "backward" + '_' +str(i)] = str(3*hsp.env_start) +':'+ str(3*hsp.env_end) 
                    except:
                        continue
            print(hmm_dict)
            file = open(path +'/'+path + '_hits.csv', 'w')
            for key, value in sorted(hmm_dict.items()):
                    file.write(str(key) +'\t' + str(value) + '\n')
                                    
def HMMread(path, record=None):
    hmm_dict = {}
    i = 0
    for infile in glob.glob(path + '/*/*.fasta'):
        try:
            qresult = SearchIO.read(infile, 'hmmer3-text')
            for i in range(len(qresult.hsps)):
                try:
                    hsp = qresult[0][i]
                    print ("HSPs")
                    print (hsp.env_start)
                    print (hsp.env_end)

                    # If we have a genome record, it means we want to pull to the start and stop codons
                    if record:
                        print ("Pull to start and stop codons")
                        start_position = getStartPostion(record, hsp.env_start)

                    hmm_dict[infile + "_" + str(i)] = str(3*hsp.env_start)+ ':' + str(3*hsp.env_end)
                    i += 1
                except:
                    continue
        except:
            continue
    return(hmm_dict)
    
#test = HMMread('hmm_outputs/DDQSY/a2')
#print(test)




