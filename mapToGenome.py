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
    queryString +=" AND genome[title] NOT shotgun[title] NOT segment[title]"
    # queryString += " AND genome[title]"

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
    genome_records = Entrez.efetch(db="nucleotide", id=queryString, rettype="gb")

    records = SeqIO.parse(genome_records, "gb")

    found_species = []
    unmappable_species = []

    for record in records:
        print (record.annotations.get('organism'))
        found_species.append(record.annotations.get('organism'))

        # print('seq dict here is ')
        # print (seqDict)

        if record.description in seqDict:

            # Prefer RefSeq sequences
            if 'RefSeq' not in record.annotations.get('keywords'):
                seqDict[record.description] = record

        else:
            seqDict[record.description] = record

    # Join all the found species together so we can quickly search to see if we didn't find something
    combined_species = '\t'.join(found_species)


    # If we didn't find a species add it to the unmappable species list
    for species in species_names:
        if species not in combined_species:
            unmappable_species.append(species)

    return unmappable_species
