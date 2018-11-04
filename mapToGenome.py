from Bio import Entrez, SeqIO
import utilities
import re
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
import gzip
from flask import flash
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import models
from http.client import IncompleteRead

Entrez.email = "gabriel.foley@uqconnect.edu.au"

# seqDict = {}
missing = []


def getSpeciesNames(seq_records, type):
    # records = SeqIO.parse(file, "fasta")
    species_names = set()
    record_ids = []

    # Make a list of the records we're querying
    for record in seq_records.values():
        # print (record.id)
        record_ids.append(record.id)
        query_string = utilities.makeQueryString(record_ids, link="+OR+")

    # Get the organism names
    protein_handle = Entrez.efetch(db=type, id=query_string, rettype="gb")
    protein_records = SeqIO.parse(protein_handle, "gb")

    for record in protein_records:
        species_names.add(record.annotations.get('organism'))

    return species_names

def getGenomeIDs(species_names):
    genome_ids = set()

    query_string = utilities.makeQueryString(species_names, "[orgn]", " OR ")

    #TODO: Let the user decide how to filter the GenBank records (i.e. let user decide if we're not accepting shotgun, segments, etc...)
    query_string_with_filters = query_string + " AND genome[title] NOT shotgun[title] NOT contig[title]  NOT segment[title] NOT plasmid[title] NOT megaplasmid[title] AND refseq[filter]"
    # query_string_with_filters = query_string + " AND genome[title] NOT shotgun[title] NOT contig[title]  NOT segment[title] NOT plasmid[title] NOT megaplasmid[title]"


    # queryString += " AND genome[title] AND project[title]"

    print('\n ***** Searching NCBI with the following query (first time)')
    print(query_string_with_filters + "\n")

    # Get a list of the genome IDs
    try:
        genome_handle = Entrez.esearch(db="genome", term= query_string, rettype="gb", idtype='acc', retmax='10000')
        # genome_handle = Entrez.esearch(db="nuccore", term= query_string_with_filters, rettype="gb", idtype='acc', retmax='10000')

        record = Entrez.read(genome_handle)

        print (record)
        print ("That was it, folks")

        for recordID in record['IdList']:
            genome_ids.add(recordID)
    except (URLError):
        pass

    return genome_ids


def get_full_genome(genome_ids):
    """
    Given an uploaded file, return the most appropriate genome record and add it to the database
    :param region_file: The file containing the regions we are interested in
    :return:
    """

    found_ids = []
    correct_alphabet_ids = []
    non_ref_seq_ids = []
    seqDict = {}

    for genome_id in genome_ids:

        genome_query_string = utilities.makeQueryString(genome_id, link="+OR+")

        print("\nThese are the genome records we identified")
        print(genome_id)
        # print (genome_query_string)
        # print (" ".join(genome for genome in genome_id) + "\n")

        check = models.GenomeRecords.query.filter_by(name=genome_id).first()

        if check :

            print ('here is the species')

            species = " ".join(models.GenomeRecords.query.filter_by(name=genome_id).first().species.split(" ")[0:2])

            print ("\n This genome record is already in the database %s" % (genome_id))
            seqDict[species] = "in_database"
        else:

            # try:
            genome_records = Entrez.elink(db="nuccore", dbfrom="genome", id=genome_id)

            print (genome_records.read())

            # genome_records = Entrez.esearch(db="genome", term=genome_id, rettype="gb")

            records = SeqIO.parse(genome_records, "gb")

            print ('here is record')

            print (records)

            found_species = []

            for record in records:
                print (record)
                print ('found something')
                found_ids.append(record.id)

                found_species.append(record.annotations.get('organism'))

                # Check if the genome is just "N" characters. Slice from the middle to account for genomes that
                # begin or end with "N" characters
                if any(nucleotide in 'A G C T' for nucleotide in record.seq[int(len(record.seq)/2):int(len(record.seq)/2 + 1000)]):
                    correct_alphabet_ids.append(record.id)

                    if record.description in seqDict:

                        # Prefer RefSeq sequences
                        if 'RefSeq' not in record.annotations.get('keywords'):
                            continue

                    else:

                        non_ref_seq_ids.append(record.id)
                        seqDict[record.annotations.get('organism')] = record

                        # seqDict[" ".join(record.annotations.get('organism').split(" ")[0:2])] = record

                # Join all the found species together so we can quickly search to see if we didn't find something
                combined_species = '\t'.join(found_species)


                # return seqDict

            # except (IncompleteRead, HTTPError, URLError):
            #     print ("We are missing genome id - ", genome_id)
            #     missing.append(genome_id)

                # pass
                # for genome_id in genome_ids:
                #     print (ex)
                #     print("Couldn't find an appropriate genome record for %s" % (genome_id))
                #     flash("Couldn't find an appropriate genome record for %s" % (genome_id))

    # Check if there were any genome IDs we couldn't find
    for genome_id in genome_ids:
        if genome_id not in found_ids:
            print("\nCouldn't find an appropriate genome record for %s" % (genome_id))
            flash("Couldn't find an appropriate genome record for %s" % (genome_id))
        if genome_id not in correct_alphabet_ids:
            print("\nThis genome record had all N characters - %s" % (genome_id))
            # flash("This genome record had all N characters - %s" % (genome_id))
        elif genome_id not in non_ref_seq_ids:
            print("\nThis genome record was omitted because another RefSeq genome for the species exists - %s" % (
            genome_id))
            flash("This genome record was omitted because another RefSeq genome for the species exists - %s" % (
            genome_id))

    return seqDict


def update_shotgun_id_dict(species_names, shotgun_id_dict):
    query_string = utilities.makeQueryString(species_names, "[orgn]", " OR ")
    query_string_with_filters = query_string + " AND genome[title] AND project[title] NOT contig[title]  NOT segment[title] NOT plasmid[title] NOT megaplasmid[title]"
    genome_ids = set()

    querypath = "tmp/tempfile"

    print ('Searching NCBI with the following query')
    print (query_string_with_filters + "\n")

    # Get a list of the genome IDs
    # try:
    genome_handle = Entrez.esearch(db="nucleotide", term= query_string_with_filters, rettype="gb", idtype='acc', retmax='10000')
    record = Entrez.read(genome_handle)
    records = ""

    for recordID in record['IdList']:
        # print (recordID)

        genome_ids.add(recordID)

        genome_query_string = utilities.makeQueryString(genome_ids, link="+OR+")

        genome_records = Entrez.efetch(db="nuccore", id=genome_query_string, rettype="gb")

        records = SeqIO.parse(genome_records, "gb")

        if records:
            # print ('found a record')
            for record in records:

                species = " ".join(record.annotations.get('organism').split()[0:2])
                strain = record.annotations['source']

                comment = record.annotations.get('comment')

                if comment:
                    # print (comment)

                    id_string = re.search('accession([\w\W].*)\.', str(comment))
                    # id_string = re.search('accession (.*)\.', comment)

                    # id_string.group(0).split("\n")
                    # print(id_string.group(0))

                    if (id_string):
                        # print (id_string.group(0).split(" ")[1])

                        # if (id_string.group(1)):

                        genome_id = id_string.group(1)

                        versionString = re.search('project[\w\W]\((.*)\)', comment)
                        version = versionString.group(1)

                        if "_" in genome_id:
                            genome_id = genome_id.split("_")[1]



                        shotgun_id_dict[genome_id.strip()] = {"version" : version, "species" : species, "source" : strain}
                        return shotgun_id_dict

                    else:
                        print("Couldn't add an appropriate record from this genome - %s" % (species_names))
                        return

        else:
            print("Couldn't add a record from this genome - %s" % (species_names))
    # except(URLError):
    #     print ("Here is a new error")
    #     return


def get_shotgun_id_dict_from_id(genome_record):
    shotgun_id_dict = {}
    genome_id = genome_record.split(".")[0]
    if "_" in genome_id:
        genome_id = genome_id.split("_")[1]
    version = genome_record.split(".")[1]

    if len(version) == 1:
        version = "0" + version

    shotgun_id_dict[genome_id] = {"version" : version}
    return shotgun_id_dict


def get_shotgun_genome(shotgun_id_dict):
    seqDict = {}

    try:
        for genome_id, genome_info in shotgun_id_dict.items():

            # print (genome_id, genome_info["version"])

            # Check if we already have this shotgun sequence

            check = models.GenomeRecords.query.filter_by(name=genome_id + " Shotgun Sequence").first()

            if check:
                print('The shotgun sequence data from %s already exists in the database' % (genome_id))

            else:

                querypath = "tmp/" + genome_id + "_shotgun_records"

                # print (genome_id)

                ftpUrl = "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/" + genome_id[0:2] + "/" + genome_id[2:4] + "/" + genome_id[0:4] + genome_info["version"] + "/" + genome_id[0:4] + genome_info["version"] + ".1.gbff.gz"
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
                print('The shotgun sequence data from %s has been written to %s \n' % (genome_id, querypath))
                new_records = SeqIO.parse(querypath, "gb")


                # Collate all of the nucleotide records together to make the genome
                shotgun_genome = ""

                for record in new_records:
                    shotgun_genome += str(record.seq)

                # Check if we have species and strain information

                species = get_shotgun_annotations(shotgun_id_dict[genome_id], "species")
                source = get_shotgun_annotations((shotgun_id_dict[genome_id]), "source")

                shotgun_seq = SeqRecord(Seq(shotgun_genome), id=genome_id + " Shotgun Sequence",
                                        annotations={"organism": species, "source": source})
                seqDict[shotgun_seq.id] = shotgun_seq

    except EOFError as error:
        flash("There was an error reading the shotgun sequence file. Please try again")
        print("There was an error reading the shotgun sequence file. Please try again")
    except URLError:
        flash("There was an error trying to read  the shotgun sequence from the FTP site. Check that it is a valid shotgun sequence identifier and try again")
        print("There was an error trying to read  the shotgun sequence from the FTP site. Check that it is a valid shotgun sequence identifier and try again")

    return seqDict


def get_shotgun_annotations(shotgun_seq, term):
    """
    Retrive annotations from within a specific shotgun sequence entry
    :param shotgun_seq: The shotgun sequence to retrieve from
    :param term: The annotation to retrieve
    :return: The value of the annotation from the shotgun sequence
    """
    if term in shotgun_seq:
        annotation = shotgun_seq[term]
    else:
        annotation = ""
    return annotation