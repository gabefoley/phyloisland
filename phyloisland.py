import requests
import lxml.etree as ET

from Bio import Entrez, SeqIO, GenBank, AlignIO, pairwise2
from Bio.SeqFeature import ExactPosition
from Bio.SubsMat.MatrixInfo import blosum62

seqDict = {}
feature_to_record = {}
selected_records = {}
selected_features = set()
unmappable = []
referenceSeqs = {}


class SmallRecord():

    def __init__(self, refseqID, emblID, uniprotID, species, sequence):
        self.refseqID = refseqID
        self.emblID = emblID
        self.uniprotID = uniprotID
        self.species = species
        self.sequence = sequence

    def __str__(self):
        return "RefSeq ID = " + self.refseqID + " EMBL ID = " + self.emblID + " UniProt ID = " + self.uniprotID

    def get_refseqID(self):
        return self.refseqID

    def get_emblID(self):
        return self.emblID

    def get_uniprotID(self):
        return self.uniprotID

    def get_species(self):
        return self.species

    def get_sequence(self):
        return self.sequence

def makeQueryString(iter, info = "", link = "", final = ""):
    queryString = ""
    for item in iter:
        queryString += item + info + link

    # Remove the final joining string from the queryString
    queryString = queryString[:-len(link)] + final

    print (queryString)
    return queryString


def getFullGenome(region_file, region_name):
    region = SeqIO.parse(region_file, "fasta")
    Entrez.email = "gabriel.foley@uqconnect.edu.au"
    species_names = set()
    genome_ids = set()

    queryString = makeQueryString(region, "+OR+")
    # queryString = ""
    #
    # for protein in region:
    #     queryString += protein.id + "+OR+"
    #
    # # Remove the final "+OR+" from the queryString
    # queryString = queryString[:-4]
    print (queryString)

    protein_handle = Entrez.efetch(db="protein", id=queryString, rettype="gb")

    protein_records = SeqIO.parse(protein_handle, "gb")

    for record in protein_records:
        species_names.add(record.annotations.get('organism'))

        #    genome_records = SeqIO.parse(genome_handle, "gb")
        #
        #     for record in genome_records:
        #         print (record)

    queryString = makeQueryString(species_names, "[orgn]", " OR ")
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

    for record in records:

        if record.description in seqDict:
            if 'RefSeq' not in record.annotations.get('keywords'):
                seqDict[record.description] = record

        # print (record.description)
        # print (record)
        # print (record.name)
        # print (record.annotations.get('keywords'))

        else:
            seqDict[record.description] = record

    print (len(seqDict))

        # print (genome_handle.read())

    # genome_records = SeqIO.parse(genome_handle, "gb")
    #
    # for record in genome_records:
    #     print(record)
    #
    # for seq_record in region:
    #     print (seq_record.id)
    #
    #     # For each protein ID get the genome ID
    #     handle = Entrez.elink(dbfrom="protein", db="genome", id =seq_record.id)
    #
    #     # print (handle.read())
    #
    #     # Store the mapping from protein ID to genome
    #     root = ET.parse(handle)
    #     for result in root.xpath('//Link/Id/text()'):
    #         prot_To_Genome[seq_record.id] = result
    #
    #     genome_handle = Entrez.esearch(db='nucleotide',
    #                                     term="Xenorhabdus bovienii[orgn] ", rettype="gb")
    #
    #     genome_records = SeqIO.parse(genome_handle, "gb")
    #
    #     for record in genome_records:
    #         print (record)

        # print (handle.read())
        # print (prot_To_Genome)

        # Get the full genome records
        # for prot, genome in prot_To_Genome.items():
        #     print (genome)
        #     genome_handle = Entrez.efetch(db="genome", id=genome, rettype="gb")
        #
        #     genome_records = SeqIO.parse(genome_handle, "gb")
        #
        #     for record in genome_records:
        #         print (record)



    # handle = Entrez.elink(dbfrom="protein", id="WP_046335775.1", db="taxonomy")
        # handle = Entrez.efetch(db="protein", id=seq_record.id, rettype="gb")
        #
        # records = SeqIO.parse(handle, "gb")
        #
        # for record in records:
        #     handle = Entrez.esearch(db='nucleotide',
        #                             term=record.annotations.get('organism') + "AND complete genome[title]")
        #     print(handle.read())




            #
        # # handle = Entrez.efetch(db="protein", id="WP_046335775.1", rettype="gb")
        #
        # records = SeqIO.parse(handle, "xml")
        #
        # for record in records:
        #     print (record)
        #     print (record.annotations.get('organism'))
        # records = SeqIO.parse("gb.gb", "gb")
        # SeqIO.write(records, "gb.gb", "gb")

        # feature_parser = GenBank.FeatureParser()
        #
        # for record in records:
        #     if small_record.get_emblID() in seqDict:
        #         seqDict[small_record.get_emblID()].annotations[set_id] = small_record.get_sequence()
        #         seqDict[small_record.get_emblID()].annotations[set_id + "_location"] = getCoords(
        #             seqDict[small_record.get_emblID()], set_id)
        #         # user_annotations.add(set_id) # User annotations now contains this ID
        #
        #     else:
        #         record.annotations[set_id] = small_record.get_sequence()
        #         record.annotations[set_id + "_location"] = getCoords(record, set_id)
        #         record.annotations[
        #             'Uniprot'] = small_record.get_uniprotID()  # Append an annotation saying this record has been identified from this set
        #         seqDict[small_record.get_emblID()] = record
        #         # user_annotations.add(set_id) # User annotations now contains this ID
        #
        #     if set_id in feature_to_record:
        #         feature_to_record[set_id].append(small_record.get_emblID())
        #     else:
        #         feature_to_record[set_id] = [small_record.get_emblID()]


def defaultValue(region1, region1_name):
        region = SeqIO.parse(region1, "fasta")

        genes = mapToGene(region, region1_name)
        # region_names.append(region1_name)

        # return genes

        for protein_id, small_record in genes.items():
            if (small_record):  # If the smallRecord mapped successfully
                getGene(protein_id, small_record, region1_name)
                getSpecies(protein_id, small_record)

            else:
                unmappable.append(protein_id)

        return "finish"


# def defaultValues(region1, region1_name):
#     region = SeqIO.parse(region1, "fasta")
#     genes = mapToGene(region, region1_name)
#     # region_names.append(region1_name)
#
#     for protein_id, small_record in genes.items():
#         if (small_record): # If the smallRecord mapped successfully
#             getGene(protein_id, small_record, region1_name)
#             getSpecies(protein_id, small_record)
#
#             # Get species
#             split_line = seqDict[small_record.get_emblID()].annotations['organism'].split()
#             species_name = split_line[0] + " " + split_line[1]
#             seqDict[small_record.get_emblID()].annotations["species"] = species_name
#         else:
#             unmappable.append(protein_id)



def getSpecies(protein_id, small_record):
    split_line = seqDict[small_record.get_emblID()].annotations['organism'].split()
    species_name = split_line[0] + " " + split_line[1]
    seqDict[small_record.get_emblID()].annotations["species"] = species_name


def mapToGene(seq_records, setID):
    protDict = {}
    # Build up a string to store multiple queries to UniProt
    queryString = ""

    for seq_record in seq_records:
        queryString += seq_record.id + "+OR+"
        protDict[seq_record.id] = "" # Add all of the protein IDs to a dictionary

    # Remove the final "+OR+" from the queryString
    queryString = queryString[:-4]

    # print (queryString)

    # Build the parameters to query UniProt with to map from RefSeq ID to EMBL database ID
    # params = {"query": queryString, "format": "tab", "columns": "database(RefSeq),database(EMBL),entry+name,protein+names,genes,database(ExpressionAtlas),organism,sequence"}
    params = {"query": queryString, "format": "tab", "columns": "database(RefSeq),database(EMBL),entry+name,organism,sequence"}

    response = requests.get("http://www.uniprot.org/uniprot/", params)

    for line in response.text.split("\n")[1:-1]: # Ignore the header line
        split_line = line.split("\t")
        print (split_line)
        if len(split_line) > 1:
            refseq_id = split_line[0][:-1] # Remove the semi-colon
            embl_id = split_line[1][:-1] # Remove the semi-colon

            uniprot_id = split_line[2] # Remove the semi-colon
            species = split_line[3]
            sequence = split_line[4]

            small_record = SmallRecord(refseq_id, embl_id, uniprot_id, species, sequence)
            # print (small_record)

            protDict[refseq_id] = small_record # Add a genbank ID for every protein ID that mapped successfully

    return protDict


def getGene(protein_id, small_record, set_id):
    Entrez.email = "gabriel.foley@uqconnect.edu.au"
    # handle = Entrez.elink(dbfrom="protein", id="WP_046335775.1", db="taxonomy")
    handle = Entrez.efetch(db="nucleotide", id= small_record.get_emblID(), rettype="gb")

    # handle = Entrez.efetch(db="protein", id="WP_046335775.1", rettype="gb")

    records = SeqIO.parse(handle, "gb")
    # records = SeqIO.parse("gb.gb", "gb")
    # SeqIO.write(records, "gb.gb", "gb")

    feature_parser = GenBank.FeatureParser()

    for record in records:
        if small_record.get_emblID() in seqDict:
            seqDict[small_record.get_emblID()].annotations[set_id] = small_record.get_sequence()
            seqDict[small_record.get_emblID()].annotations[set_id + "_location"] = getCoords(seqDict[small_record.get_emblID()], set_id)
            # user_annotations.add(set_id) # User annotations now contains this ID

        else:
            record.annotations[set_id] = small_record.get_sequence()
            record.annotations[set_id + "_location"] = getCoords(record, set_id)
            record.annotations['Uniprot'] = small_record.get_uniprotID()# Append an annotation saying this record has been identified from this set
            seqDict[small_record.get_emblID()] = record
            # user_annotations.add(set_id) # User annotations now contains this ID

        if set_id in feature_to_record:
            feature_to_record[set_id].append(small_record.get_emblID())
        else:
            feature_to_record[set_id] = [small_record.get_emblID()]


def getCoords(record, set_id):
    for feature in record.features:
        if 'translation' in feature.qualifiers:
            if feature.qualifiers['translation'][0] == record.annotations[set_id]:
                location_1_start = feature.location.start
                location_1_end = feature.location.end
                return str(location_1_start) + ":" + str(location_1_end)
    return "No coords found"


def getDistance(overlap_1, overlap_2):
    loc_1_start, loc_1_end = overlap_1.split(":")
    loc_2_start, loc_2_end = overlap_1.split(":")
    if loc_1_end > loc_2_start or loc_2_end > loc_1_start:
        return ""
    elif loc_1_end < loc_2_start:
        return str(loc_2_start - loc_1_end)
    elif loc_2_end < loc_1_start:
        return str(loc_1_start - loc_2_end)


def addToRecord(old, new, region):
    setattr(old, region, 21)
    print ("**OLD**", old)
    print ("**NEW**", new)
    print ("**REGION**", region)
    old.a2 = new.annotations[region] if region in old.annotations.keys() else "Not tested"
    old.a2_loc = new.annotations[region + "_location"] if region + "_location" in new.annotations.keys() else "Not tested"
    return old


def check_overlap(overlap_1, overlap_2):
    for record in seqDict:
        # Only check records that have both of the regions
        if overlap_1 in seqDict.get(record).annotations.keys() and overlap_2 in seqDict.get(record).annotations.keys() :
            location1 = ExactPosition(0)
            location2 = ExactPosition(0)
            for feature in seqDict.get(record).features:
                if 'translation' in feature.qualifiers:

                    # print (feature.qualifiers['translation'][0])
                    if feature.qualifiers['translation'][0] == seqDict.get(record).annotations[overlap_1]:
                        location_1_start = feature.location.start
                        location_1_end = feature.location.end
                    if feature.qualifiers['translation'][0] == seqDict.get(record).annotations[overlap_2]:
                        location_2_start = feature.location.start
                        location_2_end = feature.location.end
            if location_1_end > location_2_start or location_2_end > location_1_start:
                seqDict[record].annotations[overlap_1 + ":" + overlap_2 + " Overlap"] = "True"
                # user_annotations.add(overlap_1 + ":" + overlap_2 + " Overlap")
            elif location_1_end < location_2_start:
                seqDict[record].annotations[overlap_1 + ":" + overlap_2 + " Overlap"] = "False"
                seqDict[record].annotations[overlap_1 + ":" + overlap_2 + " Distance"] = location_2_start - location_1_end
                # user_annotations.add(overlap_1 + ":" + overlap_2 + " Overlap")
                # user_annotations.add(overlap_1 + ":" + overlap_2 + " Distance")

            elif location_2_end < location_1_start:
                seqDict[record].annotations[overlap_2 + ":" + overlap_1 + " Overlap"] = "False"
                seqDict[record].annotations[overlap_1 + ":" + overlap_2 + " Distance"] = location_1_start - location_2_end
                # user_annotations.add(overlap_2 + ":" + overlap_1 + " Overlap")
                # user_annotations.add(overlap_2 + ":" + overlap_1 + " Distance")


def checkForFeature(seqRecords, featureText):
    for record in seqRecords:
        for feature in seqRecords.get(record).features:
            if 'gene' in feature.qualifiers and 'translation' in feature.qualifiers:
                if featureText in feature.qualifiers['gene'][0]:
                    seqRecords[record].annotations[featureText] = 'Found'

                    # print ("Found", feature.qualifiers['gene'][0], feature.qualifiers['translation'][0])


def unique_strains():
    strains = set()
    unique_strains = {}
    global selected_records, referenceSeqs

    if len(referenceSeqs) < 1:
        print ("*** Warning: You haven't selected a reference sequence yet. The unique record chosen for each strain will be randomly selected ")
    else:
        for k,v in referenceSeqs.items():
            print ("Using ", k, "as the reference sequence for the", v, "feature. The unique record chosen for each strain without this feature will be randomly selected")


    for record in seqDict:

        if seqDict[record].annotations['organism'] in strains:
            pass
        else:
            strains.add(seqDict[record].annotations['organism'])
            unique_strains[record] = seqDict[record]
    selected_records = unique_strains

    # profile_menu()


def getUniqueSpecies(refseq, records):

    for record in records:
        pass

def unique_species():
    species = {}
    unique_species = {}
    global selected_records

    if len(referenceSeqs) < 1:
        print ("*** Warning: You haven't selected a reference sequence yet. The unique record chosen for each strain will be randomly selected ")
    else:
        for k,v in referenceSeqs.items():
            print ("Using ", k, "as the reference sequence for the", v, "feature. The unique record chosen for each strain without this feature will be randomly selected")

    for record in seqDict:
        if ("A1" in seqDict[record].annotations.keys()):
            # print ("LENGTH" + (str(len(unique_species))))
            # for k,v in unique_species.items():
            #     print ("here it is ", k)

            # print ("record " + record)
            split_line = seqDict[record].annotations['organism'].split()
            species_name = split_line[0] + " " + split_line[1]
            print (seqDict[record].annotations.keys())
            # print ("A1" in seqRecords[record].annotations.keys() )
            if species_name in species.keys():
                # print(species_name + " was in species")
                high_score = species[species_name][1]
                # print ("Checking alignment score for ", record, " and  ")
                new_score = getAlignmentScore(seqDict[referenceSeqs["A1"]].annotations["A1"], seqDict[record].annotations["A1"])
                # print ("New score " + str(new_score) + " and high score " + str(high_score))
                if new_score > high_score:
                    unique_species.pop(species[species_name[0]])
                    species[species_name] = [record, new_score]
                    unique_species[record] = seqDict[record]
                    # for k, v in unique_species.items():
                    #     print("here it is inside here ", k)

            else:
                # print (species_name + " was  not in species")
                high_score = getAlignmentScore(seqDict[referenceSeqs["A1"]].annotations["A1"], seqDict[record].annotations["A1"])
                # print ("GOT HERE")
                species[species_name] = [record, high_score]
                # print ("NOW HERE")
                unique_species[record] = seqDict[record]
                # print ("WHERE HERE")
                # for k, v in unique_species.items():
                #     print("here it is from here ", k)

    # for k, v in unique_species.items():
    #     print("here it is finally  ", k)
    selected_records = unique_species

    # profile_menu()


def getAlignmentScore(seq1, seq2):
    print (seq1)
    print (seq2)

    alignment_score = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5, score_only=True)
    return alignment_score

