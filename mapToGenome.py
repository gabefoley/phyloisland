from Bio import Entrez, SeqIO, GenBank, AlignIO, pairwise2
import utilities
import re
import urllib
import requests, zipfile
import io
from zipfile import ZipFile
from urllib.request import urlopen
from gzip import GzipFile


seqDict = {}


def getFullGenome(region_file):
    records = SeqIO.parse(region_file, "fasta")
    Entrez.email = "gabriel.foley@uqconnect.edu.au"
    species_names = set()
    genome_ids = set()
    record_ids = []
    # print ('returning')
    #
    # handle = Entrez.efetch(db='nucleotide', rettype='gbwithparts', id='MWTM00000000', retmode='text')
    # print (handle)
    # protein_records = SeqIO.parse(handle, "gb")
    #
    # for record in protein_records:
    #     print (record)
    # print ('*******************')
    #
    #


    # Make a list of the records we're querying
    for record in records:
        record_ids.append(record.id)
    queryString = utilities.makeQueryString(record_ids, link="+OR+")


    # Get the organism names
    #TODO: Handle not just protein sequences here

    protein_handle = Entrez.efetch(db="protein", id=queryString, rettype="gb")
    protein_records = SeqIO.parse(protein_handle, "gb")

    for record in protein_records:
        species_names.add(record.annotations.get('organism'))


    queryString = utilities.makeQueryString(species_names, "[orgn]", " OR ")


    #TODO: Let the user decide how to filter the GenBank records (i.e. let user decide if we're not accepting shotgun, segments, etc...)
    # queryString +=" AND genome[title] AND project[title] NOT contig[title]  NOT segment[title] NOT plasmid[title]"
    queryString +=" AND genome[title] NOT shotgun[title] NOT contig[title]  NOT segment[title] NOT plasmid[title]"

    # queryString += " AND genome[title] AND project[title]"

    print ('query string is')
    print (queryString)

    # Get a list of the genome IDs
    genome_handle = Entrez.esearch(db="nucleotide", term= queryString, rettype="gb", idtype='acc')
    record = Entrez.read(genome_handle)
    print (record)
    print (record['IdList'])

    for recordID in record['IdList']:
        genome_ids.add(recordID)

    # Convert GI numbers to accession numbers



    # root = ET.parse(genome_handle)
    # print (root)
    # for result in root.xpath('///Id/text()'):
    #     genome_ids.add(result)

    # for result in root.xpath('///wgs/text()'):
    #     print('result')
    #     print (result)
    queryString = utilities.makeQueryString(genome_ids, link="+OR+")
    #
    # print (queryString)
    #
    # genome_handle = Entrez.efetch(db="nuccore", term= queryString, rettype="acc")
    #
    # print (genome_handle)



    print ('this is the query string now')
    print (queryString)
    genome_records = Entrez.efetch(db="nuccore", id=queryString, rettype="gb")

    records = SeqIO.parse(genome_records, "gb")

    found_species = []
    unmappable_species = []
    idDict = {}

    for record in records:
        # print (record)
        # comment = record.annotations.get('comment')
        #
        # comment.split()
        # idString = re.search('accession (.*)\.', comment)
        #
        # id = idString.group(1)
        # versionString = re.search('project \((.*)\)', comment)
        # version = versionString.group(1)
        #
        # print ('id')
        # print (id)
        # print ('version')
        # print (version)
        #
        # if "_" in id:
        #     id = id.split("_")[1]
        #
        # print (id)
        #
        # idDict[id] = version
        #
        # for id, version in idDict.items():
        #
        #     ftpUrl = "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/" + id[0:2] + "/" + id[2:4] + "/" + id[0:4] + version + "/" + id[0:4] + version + ".1.gbff.gz"
        #     print ("********************")
        #     print (ftpUrl)
        #
        #     # Download the file from the URL
        #     zipresp = urlopen(ftpUrl)
        #     # Create a new file on the hard drive
        #     tempzip = open("tmp/tempfile.gz", "wb")
        #     # Write the contents of the downloaded file into the new file
        #     tempzip.write(zipresp.read())
        #     # Close the newly-created file
        #     tempzip.close()
        #     # Re-open the newly-created file with ZipFile()
        #     gzf = GzipFile("tmp/tempfile.gz")
        #     # Extract its contents into /tmp/mystuff
        #     # note that extractall will automatically create the path
        #     with gzf.open("tmp/tempfile.gz", 'rb') as f:
        #         file_content = f.read()
        #         print (file_content)
        #
        #
        #     # response = requests.get(ftpUrl)
        #     # zipDocument = zipfile.ZipFile(io.BytesIO(response.content))
        #     # ftpRecord = zipDocument.extractall()
        #
        #     # mysock = urllib.urlopen(ftpUrl)
        #     # memfile = io.BytesIO(mysock.read())
        #     # with ZipFile(memfile, 'r') as myzip:
        #     #     f = myzip.open()
        #     #     content = f.read()
        #     #     print (content)
        #
        #     # ftpRecord = urllib.request.urlopen(ftpUrl)
        #     # print(ftpRecord.read(100))
        #
        # print (ftpUrl)
        print (record.annotations.get('organism'))
        found_species.append(record.annotations.get('organism'))
        print (record.seq.alphabet)
        print (record.seq[0:1000])
        if any(nucleotide in 'A G C T' for nucleotide in record.seq[0:10000]):

            if record.description in seqDict:

                # Prefer RefSeq sequences
                if 'RefSeq' not in record.annotations.get('keywords'):
                    continue
                    # seqDict[record.description] = record

            else:
                seqDict[record.description] = record

        # Join all the found species together so we can quickly search to see if we didn't find something
        combined_species = '\t'.join(found_species)


        # If we didn't find a species add it to the unmappable species list
        for species in species_names:
            if species not in combined_species:
                unmappable_species.append(species)

    return unmappable_species
