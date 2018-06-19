# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 09:39:05 2018

@author: Rowan
"""

import os
import glob
from Bio import SearchIO

paths = ['prof1', 'tcda1', 'tcda2','tcda4','xpta1','xpta2','yena1','yena2']


max_length = 0
hmm_dict = []
for path in paths:
    for infile in glob.glob( os.path.join(path, '*.fasta')):
        qresult = SearchIO.read(infile, 'hmmer3-text')
        maxscore = 0
        for i in range(len(qresult.hsps)):
                try:
                    hsp = qresult[0][i]
                    hit = qresult[0]
                    hmm_dict.append([str(hsp.env_start), str(hsp.env_end), str(hsp.bitscore), infile])
                except:
                    continue
            

print(hmm_dict)
"""
for hit in hmm_dict.keys():
    if hmm_dict[hit] != []:
        new_hmm_dict[hit] = hmm_dict[hit]
        

for hit in new_hmm_dict.keys():
    if len(new_hmm_dict[hit]) < max_length:
        while len(new_hmm_dict[hit]) < max_length:
            new_hmm_dict[hit].append(['-', '-', '-', '-', '-'])
    #if 
"""







import csv
file = open('testing.csv', 'w')
with file:
    writer = csv.writer(file)
    writer.writerows(hmm_dict)