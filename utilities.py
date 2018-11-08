import os
import glob
from flask import flash
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet


def makeQueryString(iter, info = "", link = "", final = ""):
    """
    Take a list of items and concatanate them together
    :param iter: List of items to concatanate
    :param info: Additional piece of info to add behind each item
    :param link: Linking phrase to link all items
    :param final: Final phrase to appened to the enf of the string
    :return: String formatted with items and linking phrases
    """
    queryString = ""

    # If it is a single record passed as a string, don't cut it up
    if isinstance(iter, str):
        queryString += iter + info + link

    else:
        for item in iter:
            queryString += item + info + link

    # Remove the final joining string from the queryString
    queryString = queryString[:-len(link)] + final

    return queryString

def removeFile(*args):
    """
    Remove files in the list from the directory

    :param args: Files to remove
    :return:
    """
    for arg in args:
        os.remove(arg)


def createFASTAFromHMMOutput(seq_record, hmmerout, species, region, amino_acid=True):
    seq_list = []
    organism = seq_record.annotations['organism']
    for result, position in hmmerout[0].items():
        start = int(position.split(":")[0])
        end = int(position.split(":")[1])
        formatted_id = "_".join([species, organism, region, position]).replace(" ", "_")
        if 'forward' in result:
            if amino_acid:
                seq = SeqRecord(Seq(str(Seq.translate(seq_record.seq[start:end]))), id=formatted_id, description="")
            else:
                seq = SeqRecord(Seq(str(seq_record.seq[start:end])), id=formatted_id, description="")
        elif 'backward' in result:
            if amino_acid:
                seq = SeqRecord(Seq.translate(Seq(str(seq_record.seq[start:end])).reverse_complement()), id=formatted_id, description="")
            else:
                seq = SeqRecord(Seq(str(seq_record.seq[start:end])).reverse_complement(), id=formatted_id, description="")
        seq_list.append(seq)
    return seq_list

def saveFASTA(align_list, filepath):
    SeqIO.write(align_list, filepath, "fasta")
    print ("FASTA file was saved to %s" % (filepath))
