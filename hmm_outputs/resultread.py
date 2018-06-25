# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 09:39:05 2018

@author: Rowan
"""

import os
import glob
from Bio import SearchIO
import csv
import pandas as pd

paths = ['DDQSY', 'FJWHO', 'ITMTT']



def resultRead(paths):
    hmm_dict = {}
    for path in paths:
        for infile in glob.glob( os.path.join(path, '*.fasta')):
            qresult = SearchIO.read(infile, 'hmmer3-text')
            for i in range(len(qresult.hsps)):
                    try:
                        hsp = qresult[0][i]
                        hmm_dict[infile[6:-24] + '_'+str(i)] = str(hsp.env_start) +':'+ str(hsp.env_end)
                    except:
                        continue
            print(hmm_dict)
            file = open('testing.csv', 'w')
            for key, value in sorted(hmm_dict.items()):
                    file.write(str(key) +'\t' + str(value) + '\n')
                                    
resultRead(paths)