from Bio import Entrez, SeqIO, GenBank, AlignIO, pairwise2
import utilities
import re
import urllib
import requests, zipfile
import io
from zipfile import ZipFile
from urllib.request import urlopen
import gzip
import os
import time
from flask import flash


seqDict = {}


def getFullGenome(region_file):
    """
    Given an uploaded file, return the most appropriate genome record and add it to the database
    :param region_file: The file containing the regions we are interested in
    :return:
    """

    records = SeqIO.parse(region_file, "fasta")
    Entrez.email = "gabriel.foley@uqconnect.edu.au"
    species_names = set()
    genome_ids = set()
    record_ids = []


    # Make a list of the records we're querying
    for record in records:
        record_ids.append(record.id)
        query_string = utilities.makeQueryString(record_ids, link="+OR+")


    # Get the organism names
    protein_handle = Entrez.efetch(db="protein", id=query_string, rettype="gb")
    protein_records = SeqIO.parse(protein_handle, "gb")

    for record in protein_records:
        species_names.add(record.annotations.get('organism'))

    query_string = utilities.makeQueryString(species_names, "[orgn]", " OR ")


    #TODO: Let the user decide how to filter the GenBank records (i.e. let user decide if we're not accepting shotgun, segments, etc...)
    query_string_with_filters = query_string + " AND genome[title] NOT shotgun[title] NOT contig[title]  NOT segment[title] NOT plasmid[title]"

    # queryString += " AND genome[title] AND project[title]"

    print ('Searching NCBI with the following query')
    print (query_string_with_filters + "\n")

    # Get a list of the genome IDs
    genome_handle = Entrez.esearch(db="nucleotide", term= query_string_with_filters, rettype="gb", idtype='acc')
    record = Entrez.read(genome_handle)

    for recordID in record['IdList']:
        genome_ids.add(recordID)

    genome_query_string = utilities.makeQueryString(genome_ids, link="+OR+")

    if (genome_query_string == ""):
        print ("We didn't identify any genome records. Attempting to search for shotgun sequenced genomes \n")
        getShotgunGenome(query_string)


    else:
        print ("These are the genome records we identified")
        print (" ".join(genome for genome in genome_ids) + "\n")

    genome_records = Entrez.efetch(db="nuccore", id=genome_query_string, rettype="gb")

    records = SeqIO.parse(genome_records, "gb")

    found_species = []
    unmappable_species = []

    for record in records:

        found_species.append(record.annotations.get('organism'))

        # Check if the genome is just "N" characters. Slice from the middle to account for genomes that begin or
        # end with "N" characters
        if any(nucleotide in 'A G C T' for nucleotide in record.seq[int(len(record.seq)/2):int(len(record.seq)/2 + 1000)]):

            if record.description in seqDict:

                # Prefer RefSeq sequences
                if 'RefSeq' not in record.annotations.get('keywords'):
                    continue

            else:
                seqDict[record.description] = record

        # Join all the found species together so we can quickly search to see if we didn't find something
        combined_species = '\t'.join(found_species)

        # If we didn't find a species add it to the unmappable species list
        for species in species_names:
            if species not in combined_species:
                unmappable_species.append(species)

    return unmappable_species


def getShotgunGenome(query_string):
    query_string_with_filters = query_string + " AND genome[title] AND project[title] NOT contig[title]  NOT segment[title] NOT plasmid[title]"
    idDict = {}
    genome_ids = set()
    querypath = "tmp/tempfile"

    print ('Searching NCBI with the following query')
    print (query_string_with_filters + "\n")

    # Get a list of the genome IDs
    genome_handle = Entrez.esearch(db="nucleotide", term= query_string_with_filters, rettype="gb", idtype='acc')
    record = Entrez.read(genome_handle)

    for recordID in record['IdList']:
        genome_ids.add(recordID)

        genome_query_string = utilities.makeQueryString(genome_ids, link="+OR+")

        genome_records = Entrez.efetch(db="nuccore", id=genome_query_string, rettype="gb")

        records = SeqIO.parse(genome_records, "gb")

        found_species = []
        unmappable_species = []

    for record in records:
        comment = record.annotations.get('comment')

        comment.split()
        idString = re.search('accession (.*)\.', comment)

        id = idString.group(1)
        versionString = re.search('project \((.*)\)', comment)
        version = versionString.group(1)


        if "_" in id:
            id = id.split("_")[1]

        idDict[id] = version

        for id, version in idDict.items():

            ftpUrl = "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/" + id[0:2] + "/" + id[2:4] + "/" + id[0:4] + version + "/" + id[0:4] + version + ".1.gbff.gz"
            print ("********************")
            print (ftpUrl)

            # Download the file from the URL
            zipresp = urlopen(ftpUrl)
            # Create a new file on the hard drive
            tempzip = open(querypath + ".gz", "wb")
            # Write the contents of the downloaded file into the new file
            tempzip.write(zipresp.read())
            # Close the newly-created file
            tempzip.close()
            # Re-open the newly-created file

            file_from_zip = gzip.open(querypath + ".gz", mode="rb")

            with open(querypath, 'w') as query_file:
                for line in file_from_zip:
                    query_file.write(line.decode('utf-8'))

            file_from_zip.close()
            new_records = SeqIO.parse(querypath, "gb")

            # for record in new_records:
            #     print(record.id)
            #     print (record.seq)



            # with gzip.open("tmp/tempfilenew.gz", 'rb') as f:
            #     file_content = f.read().decode('utf-8')
            #     print ('fono')
            #     new_record = SeqIO.read(file_content, "gb")
            #     print ('bono')
            #     print (new_record.name)



