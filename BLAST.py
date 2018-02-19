from Bio import SeqIO
from Bio.Blast.Applications import NcbipsiblastCommandline, NcbitblastnCommandline
from Bio.Blast import NCBIXML
import sys
import subprocess
import os


def makeBlastDB(records, dbpath):
    """
    Make a local BLAST database from selected records on the command line
    :param records:
    :return:
    """
    SeqIO.write(records, dbpath, "fasta")
    blastdb_cmd = 'makeblastdb -in %s -dbtype nucl -title temp_blastdb' % (dbpath)
    DB_process = subprocess.Popen(blastdb_cmd,
                              shell=True,
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
    DB_process.wait()
    stdout, stderr = DB_process.communicate()
    return stdout

def psiBlast(dbFile, queryFile, outfile, evalNum):
    psi_cline = NcbipsiblastCommandline('psiblast', db = dbFile , query = queryFile , evalue = evalNum , out =  outfile, outfmt = 5, out_pssm = "_pssm")
    p = subprocess.Popen(str(psi_cline),stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=(sys.platform!="win32"))

def tBlastN(dbFile, queryFile, outfile, evalNum):

    tN_cline = NcbitblastnCommandline(db = dbFile , query = queryFile , evalue = evalNum , out = outfile, outfmt = 5)
    p = subprocess.Popen(str(tN_cline),stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=(sys.platform!="win32"))


def getBlastInfo(xmlFile, closest_to=-1):
    """
    Get the location and sequence of the top scoring hits from a BLAST results
    :param xmlFile: BLAST records
    :return:
    """
    blast_info = {}
    blast_parser = NCBIXML.parse(xmlFile)

    try:
        for record in blast_parser:
            for alignment in record.alignments:
                print ("BLAST search found %s hits " % (len(alignment.hsps)))

                # If we just want the top hit and don't care about proximity
                if closest_to == -1:
                    for hsp in alignment.hsps:
                        # Check if this is on the reverse strand
                        if (hsp.frame[1] < 0):
                            blast_info["location"] = str(alignment.length - hsp.sbjct_end) + ":" + \
                                                     str( alignment.length - hsp.sbjct_start)
                        else:
                            blast_info["location"] = str(hsp.sbjct_start) + ":" + str(hsp.sbjct_end)
                        blast_info["sequence"] = hsp.sbjct
                        break
                # Otherwise we have a genomic location we want to get closest to
                else:
                    # Set distance as high as possible
                    distance = alignment.length
                    for hsp in alignment.hsps:
                        # Calculate the distance between this high scoring pair and the region we want to get closest to
                        hsp_distance = hsp.sbjct_start - int(closest_to)

                        # If this distance is shorter then set this as the closest pair
                        if hsp_distance < distance:
                            distance = hsp_distance
                            # Check if this is on the reverse strands
                            if (hsp.frame[1] < 0):
                                blast_info["location"] = str(alignment.length - hsp.sbjct_end)  + ":" + \
                                                         str(alignment.length - hsp.sbjct_start)
                            else:
                                blast_info["location"] = str(hsp.sbjct_start) + ":" + str(hsp.sbjct_end)
                            blast_info["sequence"] = hsp.sbjct
    except ValueError:
        print ("Couldn't find a BLAST hit")


    return blast_info

