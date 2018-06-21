# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 09:39:05 2018

@author: Rowan
"""

import os
import glob
from Bio import SearchIO
import csv

paths = ['prof1', 'tcda1', 'tcda2','tcda4','xpta1','xpta2','yena1','yena2']



def resultRead(paths):
    hmm_dict = {}
    for path in paths:
        for infile in glob.glob( os.path.join(path, '*.fasta')):
            qresult = SearchIO.read(infile, 'hmmer3-text')
            for i in range(len(qresult.hsps)):
                    try:
                        hsp = qresult[0][i]
                        hmm_dict[str(qresult.accession)] = [str(hsp.env_start), str(hsp.env_end)]
                    except:
                        continue
            file = open('testing.csv', 'w')
            with file:
                writer = csv.writer(file)
                writer.writerows(hmm_dict)