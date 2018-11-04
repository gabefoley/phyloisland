# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 09:39:05 2018

@author: Rowan, Gabe
"""

import os
import glob
from Bio import SearchIO
from Bio.Seq import Seq
import re


def expandStartPostion(record, hit_start, strand, track_forwards = True):

    if strand == "forward":
        codon_list = ["TGA", "TAA", "TAG"]
    elif strand == "backward":
        codon_list = ["TCA", "TTA", "CTA"]

    codon = ""
    first_pos = hit_start
    second_pos = hit_start + 3

    start_pos = first_pos # variable to keep track of all the potential start positions (methionines) we run into
    while  first_pos -3 > 0:
        codon = record.sequence[first_pos: second_pos]

        if codon == "ATG":
            start_pos = first_pos

        if codon in codon_list:
            break

        second_pos = first_pos
        first_pos = first_pos - 3

    if strand == "forward":
        # If we want to also look forwards and our start position isn't a methionine
        if track_forwards and start_pos != "ATG":
            first_pos = start_pos
            second_pos = start_pos + 3

            while first_pos + 3 < len(record.sequence):
                codon = record.sequence[first_pos: second_pos]
                if codon == "ATG":
                    start_pos = first_pos
                    break

                first_pos = second_pos
                second_pos = first_pos + 3

            return start_pos


        return start_pos
    else:
        return first_pos + 3





def expandEndPosition(record, hit_end, strand, track_forwards=True):

    if strand == "forward":
        codon_list = ["TGA", "TAA", "TAG"]
    elif strand == "backward":
        codon_list = ["TCA", "TTA", "CTA"]

    codon = ""
    first_pos = hit_end - 3
    second_pos = hit_end

    start_pos = second_pos # variable to keep track of all the potential start positions (methionines) we run into

    while  first_pos +3 < len(record.sequence):
        codon = record.sequence[first_pos: second_pos]
        if codon == "CAT":
            start_pos = second_pos
        if codon in codon_list:
            break
        first_pos = second_pos
        second_pos = first_pos + 3

    if strand == "forward":
        return second_pos - 3
    else:
        # If we want to also look forwards and our start position isn't a methionine
        if track_forwards and start_pos != "CAT":
            first_pos = start_pos
            second_pos = start_pos - 3

            while first_pos - 3 > 0:
                codon = record.sequence[first_pos: second_pos]
                if codon == "CAT":
                    start_pos = first_pos
                    break

                first_pos = second_pos
                second_pos = first_pos - 3

            return start_pos
        return start_pos




def HMMread(path, record=None, expand=False):
    # print('here we are in hmmRead')
    hmm_dict = {}
    i = 0
    for infile in glob.glob(path + '/*/*.fasta'):
        try:
            qresult = SearchIO.read(infile, 'hmmer3-text')
            strand_regex = re.search(r'_.{3,4}ward_\d', infile)

            if strand_regex:
                frame = strand_regex.group()[-1]
                strand = strand_regex.group().split("_")[1]

                correction = 0
                if strand == "forward":
                    if frame == "0":
                        correction = 2
                    elif frame == "1":
                        correction = 1
                elif strand == "backward":
                    if frame == "0":
                        correction = 1
                    elif frame == "2":
                        correction = -1

            for i in range(len(qresult.hsps)):
                try:
                    hsp = qresult[0][i]



                    if strand == "forward":
                        start = ((hsp.hit_start + 1) * 3) - correction - 1
                        end = ((hsp.hit_end + 1) * 3) - correction - 1

                    # If on the backwards strand we need to update the positions to have them in correct 5' to 3'
                    elif strand == "backward":
                        start = len(record.sequence) + correction - (hsp.hit_end * 3) - 1
                        end = len(record.sequence) + correction  - (hsp.hit_start * 3) - 1

                    # If we have a genome record, it means we want to pull to the start and stop codons
                    if expand:

                        start = expandStartPostion(record, start, strand)
                        end = expandEndPosition(record, end, strand)

                    # After expanding, we might have the exact region already identified - don't add multiple regions in
                    if str(start) + ":" + str(end) in hmm_dict.values():
                        print ("Found two identical regions, skipping adding this record %s at position %s : %s" % (infile + "_" + str(i), str(start), str(end)))

                    else:
                        hmm_dict[infile + "_" + str(i)] = str(start) + ':' + str(end)
                    i += 1
                except ValueError:
                    continue
        except ValueError:
            continue
    return (hmm_dict)






