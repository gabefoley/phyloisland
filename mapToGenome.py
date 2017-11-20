from Bio import Entrez, SeqIO, GenBank, AlignIO, pairwise2
import utilities
import lxml.etree as ET


seqDict = {}


def getFullGenome(region_file, region_name):
    records = SeqIO.parse(region_file, "fasta")
    Entrez.email = "gabriel.foley@uqconnect.edu.au"
    species_names = set()
    genome_ids = set()
    record_ids = []



    for record in records:
        record_ids.append(record.id)
    queryString = utilities.makeQueryString(record_ids, link="+OR+")
    # queryString = ""
    #
    # for protein in records:
    #     queryString += protein.id + "+OR+"
    #
    # # Remove the final "+OR+" from the queryString
    # queryString = queryString[:-4]

    protein_handle = Entrez.efetch(db="protein", id=queryString, rettype="gb")

    protein_records = SeqIO.parse(protein_handle, "gb")

    for record in protein_records:
        species_names.add(record.annotations.get('organism'))

        #    genome_records = SeqIO.parse(genome_handle, "gb")
        #
        #     for record in genome_records:
        #         print (record)

    queryString = utilities.makeQueryString(species_names, "[orgn]", " OR ")
    # queryString = ""
    # for species in species_names:
    #     queryString += species + "[orgn] OR "
    #
    # # Remove the final "+OR+" from the queryString
    # queryString = queryString[:-4]

    queryString +=" AND genome[title] NOT shotgun[title] NOT segment[title]"
    print (queryString)


    genome_handle = Entrez.esearch(db="nucleotide", term= queryString, rettype="gb")

    root = ET.parse(genome_handle)
    for result in root.xpath('///Id/text()'):
        genome_ids.add(result)

    queryString = ""

    for genome_id in genome_ids:
        queryString += genome_id + "+OR+"

    genome_records = Entrez.efetch(db="nucleotide", id=queryString, rettype="gb")

    records = SeqIO.parse(genome_records, "gb")

    found_species = []
    unmappable_species = []

    for record in records:
        found_species.append(record.annotations.get('organism'))

        if record.description in seqDict:
            if 'RefSeq' not in record.annotations.get('keywords'):
                seqDict[record.description] = record

        # print (record.description)
        # print (record)
        # print (record.name)
        # print (record.annotations.get('keywords'))

        else:
            seqDict[record.description] = record

    for species in species_names:
        if species not in found_species:
            unmappable_species.append(species)

    return unmappable_species
