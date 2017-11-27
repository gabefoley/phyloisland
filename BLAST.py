from Bio import SeqIO
from Bio.Blast.Applications import NcbipsiblastCommandline, NcbitblastnCommandline
from Bio.Blast import NCBIXML

import sys
import subprocess
import os


def makeBlastDB(records):
    SeqIO.write(records, 'files/temp_blastfiles.fasta', "fasta")
    blastdb_cmd = 'makeblastdb -in files/temp_blastfiles.fasta  -dbtype nucl -title temp_blastdb'
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
    print ("got to tBlastN")
    tN_cline = NcbitblastnCommandline(db = dbFile , query = queryFile , evalue = evalNum , out =  "files/psi.xml", outfmt = 5)
    p = subprocess.Popen(str(tN_cline),stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=(sys.platform!="win32"))


def getBlastLocation(xmlFile):
    blast_parser = NCBIXML.parse(xmlFile)
    for record in blast_parser:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                return (str(hsp.sbjct_start) + ":" + str(hsp.sbjct_end))


def parseBlastResults(xmlFile):
    print ("got to pbr")
    blast_parser = NCBIXML.parse(xmlFile)
    for record in blast_parser:
        print (record)
        print (record.alignments)
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print (hsp.sbjct_start)
                print (hsp.sbjct_end)
                print ()
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('score:', hsp.score)
                print('gaps:', hsp.gaps)
                print('e-value:', hsp.expect)
                print(hsp.query[0:90] +'...')
                print(hsp.match[0:90] +'...')
                # print(hsp.subject[0:90] +'...')


def removeFile(*args):
    for arg in args:
        os.remove(arg)


# seqs = SeqIO.parse("files/YenA1_single.fasta", "fasta")
# makeBlastDB(seqs)
# queryBlastDB('files/temp_blastfiles.fasta', "files/YenA1_single.fasta", 0.00005)
# parseBlastResults(open("files/psi.xml"))
# removeFile("files/temp_blastfiles.fasta", "psi_xml")
