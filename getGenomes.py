import argparse, subprocess, sys, fnmatch
from ftplib import FTP
import gzip
from Bio import SeqIO
from datetime import datetime, MINYEAR

def get_latest_file(seq_records, folder, file_type="gff", database="refseq"):

        print("Getting single")

        print ("rsync --list-only %s rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/%s/*/ %s" % (file_type,
                                                                                                   database,
                                                                                                   seq_records.replace(
                                                                                                       " ", "_"),
                                                                                                   folder, "./tmp"))

        process = subprocess.Popen(
            "rsync --list-only %s rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/%s/*/ %s" % (file_type,
                                                                                                   database,
                                                                                                   seq_records.replace(
                                                                                                       " ", "_"),
                                                                                                   folder, "./tmp"),
            stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

        out, err = process.communicate()
        errcode = process.returncode

        print("errcode ", errcode)

        print (out.decode("utf-8"))

        latest = ""
        latest_time = datetime.strptime("1900", '%Y')

        for line in out.decode("utf-8").split("\n"):
            if line.startswith("-r--r--r--"):
                print (line)
                splitline = line.split(" ")
                splitline = [x for x in splitline if x!= ""]
                print (splitline)
                date_time = " ".join(splitline[2:4])



                datetime_object = datetime.strptime(date_time, '%Y/%m/%d %X')

                if datetime_object > latest_time:
                    latest_time = datetime_object
                    latest = splitline[4]

        if errcode != 0:
            return 'Fail'
        else:
            return latest


def add_genome(seq_records, folder, database="refseq", single=True):

    file_type = "--exclude='*cds_from*' --exclude='*rna_from*' --include='*genomic.fna.gz' --exclude='*'"

    if single:
        filename = get_latest_file(seq_records, folder, file_type=file_type, database=database)
        if filename == "Fail":
            return "Fail"
    else:
        filename = ""

    try:
        process = subprocess.Popen(
                "rsync -Lrtv --chmod=+rwx -p %s rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/%s/*/%s %s" % (
                    file_type, database, seq_records.replace(" ", "_"), folder, filename, "./tmp"), stderr=subprocess.PIPE,
                stdout=subprocess.PIPE, shell=True)

        out, err = process.communicate()
        errcode = process.returncode

        # output = process.stdout

        print ("errcode ", errcode)
        print ("output ", out.decode("utf-8"))

        if errcode != 0:
            return 'Fail'



        file_list = out.decode("utf-8").split("receiving file list ... done")[1].split("sent")[0]

        if not file_list:
            return "Fail"

        for filename in file_list.split("\n"):
            if len(filename) > 2:
                print (filename)
                filepath = "./tmp/" + filename

                file_from_zip = gzip.open(filepath, mode="rb")

                outpath = ".".join(filepath.split(".")[0:-1]) + ".fasta"
                print (outpath)

                print (file_from_zip)

                # with open(outpath, 'w') as query_file:
                #     for line in file_from_zip:
                #         while "plasmid" not in line.decode('utf-8'):
                #             query_file.write(line.decode('utf-8'))
                #
                # file_from_zip.close()
                # print('The shotgun sequence data has been written to %s \n' % (outpath))
                # new_records = SeqIO.parse(filepath, "gb")
            return

    except subprocess.CalledProcessError as exc:
        return "Fail"


