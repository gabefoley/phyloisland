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

def psiBlast(dbFile, queryFile, evalNum):
    psi_cline = NcbipsiblastCommandline('psiblast', db = dbFile , query = queryFile , evalue = evalNum , out =  "files/psi.xml", outfmt = 5, out_pssm = "_pssm")
    p = subprocess.Popen(str(psi_cline),stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=(sys.platform!="win32"))

def tBlastN(dbFile, queryFile, evalNum):
    tN_cline = NcbitblastnCommandline(db = dbFile , query = queryFile , evalue = evalNum , out =  "files/psi.xml", outfmt = 5)
    p = subprocess.Popen(str(tN_cline),stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=(sys.platform!="win32"))


def getBlastLocation(xmlFile):
    """
    Get the location of the top scoring hits from a BLAST results
    :param xmlFile:
    :return:
    """
    blast_parser = NCBIXML.parse(xmlFile)
    for record in blast_parser:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                return (str(hsp.sbjct_start) + ":" + str(hsp.sbjct_end))