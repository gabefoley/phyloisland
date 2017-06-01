# import sys
import requests
import subprocess
import re
# from collections import defaultdict
# from io import StringIO

from Bio import Entrez, SeqIO, GenBank, AlignIO, pairwise2
from Bio.SeqFeature import ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.SubsMat.MatrixInfo import blosum62
from subprocess import DEVNULL






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

def defaultValue(region1, region1_name):
        region = SeqIO.parse(region1, "fasta")


        genes = mapToGene(region, region1_name)
        region_names.append(region1_name)

        # return genes

        for protein_id, small_record in genes.items():
            if (small_record):  # If the smallRecord mapped successfully
                getGene(protein_id, small_record, region1_name)
                getSpecies(protein_id, small_record)

            else:
                unmappable.append(protein_id)

        return "finish"


def defaultValues(region1, region1_name):
    region = SeqIO.parse(region1, "fasta")

    genes = mapToGene(region, region1_name)
    region_names.append(region1_name)



    for protein_id, small_record in genes.items():
        if (small_record): # If the smallRecord mapped successfully
            getGene(protein_id, small_record, region1_name)
            getSpecies(protein_id, small_record)

            # Get species
            split_line = seqDict[small_record.get_emblID()].annotations['organism'].split()
            species_name = split_line[0] + " " + split_line[1]
            seqDict[small_record.get_emblID()].annotations["species"] = species_name
        else:
            unmappable.append(protein_id)


    # region = SeqIO.parse(region2, "fasta")
    # genes = mapToGene(region, region2_name)
    # region_names.append(region2_name)
    #
    # for protein_id, small_record in genes.items():
    #     if (small_record):  # If the smallRecord mapped successfully
    #         getGene(protein_id, small_record, region2_name)
    #     else:
    #         unmappable.append(protein_id)


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
            user_annotations.add(set_id) # User annotations now contains this ID

        else:
            record.annotations[set_id] = small_record.get_sequence()
            record.annotations[set_id + "_location"] = getCoords(record, set_id)
            record.annotations['Uniprot'] = small_record.get_uniprotID()# Append an annotation saying this record has been identified from this set
            seqDict[small_record.get_emblID()] = record
            user_annotations.add(set_id) # User annotations now contains this ID

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
                user_annotations.add(overlap_1 + ":" + overlap_2 + " Overlap")
            elif location_1_end < location_2_start:
                seqDict[record].annotations[overlap_1 + ":" + overlap_2 + " Overlap"] = "False"
                seqDict[record].annotations[overlap_1 + ":" + overlap_2 + " Distance"] = location_2_start - location_1_end
                user_annotations.add(overlap_1 + ":" + overlap_2 + " Overlap")
                user_annotations.add(overlap_1 + ":" + overlap_2 + " Distance")

            elif location_2_end < location_1_start:
                seqDict[record].annotations[overlap_2 + ":" + overlap_1 + " Overlap"] = "False"
                seqDict[record].annotations[overlap_1 + ":" + overlap_2 + " Distance"] = location_1_start - location_2_end
                user_annotations.add(overlap_2 + ":" + overlap_1 + " Overlap")
                user_annotations.add(overlap_2 + ":" + overlap_1 + " Distance")







def checkForFeature(seqRecords, featureText):
    for record in seqRecords:
        for feature in seqRecords.get(record).features:
            if 'gene' in feature.qualifiers and 'translation' in feature.qualifiers:
                if featureText in feature.qualifiers['gene'][0]:
                    seqRecords[record].annotations[featureText] = 'Found'

                    # print ("Found", feature.qualifiers['gene'][0], feature.qualifiers['translation'][0])

def checkForFlankingFeature():
    pass

def getFullGenome(seqRecords):
    pass

def recordContains():
    pass

def getUniqueSpecies():
    pass


def print_menu():  ## Your menu design here
    print (30 * "-", "MENU", 30 * "-")
    print ("1. Add a region of interest")
    print ("2. Map a region of interest *** NOT IMPLEMENTED ***")
    print ("3. Select records and build a profile *** NOT IMPLEMENTED *** ")
    print ("4. Print summary of existing gene records")
    print ("5. Exit")
    print (67 * "-")

def print_map_menu():
    print (30 * "-", "MENU", 30 * "-")
    print ("1. Check for presence of regions")
    print ("2. Check for overlap of regions")
    print ("3. Return to main menu")
    print (67 * "-")



# loop = True
#
# while loop:  ## While loop which will keep going until loop = False
#     print_menu()  ## Displays menu
#     choice = input("Enter your choice [1-5]: ")
#
#     if choice == '1':
#         filename = input("Enter the filename: ")
#         region_name = input("Enter a name for this region: ")
#         filename = "files/YenA1_tiny.fasta"
#         region = SeqIO.parse(filename, "fasta")
#         genes = mapToGene(region, region_name)
#
#         for protein_id, small_record in genes.items():
#             if (small_record):  # If the smallRecord mapped successfully
#                 getGene(protein_id, small_record, region_name)
#             else:
#                 unmappable.append(protein_id)
#
#
#     elif choice == '2':
#         print_map_menu()
#     elif choice == '3':
#         print ("Menu 3 has been selected")
#         ## You can add your code or functions here
#     elif choice == '4':
#         for record in seqRecords:
#             print (seqRecords[record])
#         ## You can add your code or functions here
#     elif choice == '5':
#         print ("Exiting")
#         ## You can add your code or functions here
#         loop = False  # This will make the while loop to end as not value of loop is set to False
#     else:
#         # Any integer inputs other than values 1-5 we print an error message
#         input("Wrong option selection. Enter any key to try again..")

# Import the modules needed to run the script.
import sys, os

# Main definition - constants
menu_actions = {}
annotate_dict = {'1': '6', '2' : 'main_menu'}
profile_dict = {'1' : '18', '2': '10', '3' : '11', '4' : '4', '5' : '13','6': '19', '7' : 'main_menu'}
record_dict = {'1' : '14', '2' : '15', '3' : 'back'}
summary_dict = {'1' : '16', '2' : '17', '3' : 'main_menu'}

region_names = []
user_annotations = set()




# =======================
#     MENUS FUNCTIONS
# =======================

# Main menu
def main_menu():
    os.system('clear')
    print (30 * "-", "MAIN MENU", 30 * "-")
    print ("1. Map a feature to a record")
    print ("2. Annotate a record")
    print ("3. Select records and build a profile")
    print ("4. Print summary of gene records")
    print ("5. Exit")
    print (71 * "-")
    choice = input("Enter your choice [1-5]: ")
    exec_menu(choice)

    return


# Execute menu1

def exec_menu(choice):
    # os.system('clear')
    ch = choice.lower()
    if ch == '':
        menu_actions['main_menu']()
    else:
        try:
            menu_actions[ch]()
        except KeyError:
            print
            "Invalid selection, please try again.\n"
            menu_actions['main_menu']()
    return

def map_to_record():

    try:
        filename = input("Enter the filename: ")
        region = SeqIO.parse(filename, "fasta")
        region_name = input("Enter a name for this region: ")
        genes = mapToGene(region, region_name)
        region_names.append(region_name)

    except FileNotFoundError:
        print ("\n *** ERROR: Filename doesn't exist. *** \n ")
        main_menu()
        return


    for protein_id, small_record in genes.items():
        if (small_record):  # If the smallRecord mapped successfully
            getGene(protein_id, small_record, region_name)
        else:
            unmappable.append(protein_id)
    main_menu()


def annotate_menu():
    os.system('clear')
    print (30 * "-", "ANNOTATE RECORDS", 30 * "-")
    print ("1. Check for overlap of regions")
    print ("2. Return to main menu")
    print (67 * "-")
    choice = input(" >>  ")
    choice = annotate_dict[choice]
    exec_menu(choice)
    return


def overlap_menu():
    print ("Check overlap")
    overlap_1 = input("First region to check for overlap ")

    if overlap_1 not in region_names:
        print ("\n *** ERROR:", overlap_1, "is not a known region. *** \n")
        annotate_menu()
        return
    overlap_2 = input("Second region to check for overlap ")

    if overlap_2 not in region_names:
        print ("\n *** ERROR:", overlap_2, "is not a know region. *** \n")
        annotate_menu()
        return

    check_overlap(overlap_1, overlap_2)

    annotate_menu()

def profile_menu():
    os.system('clear')
    print (20 * "-", " SELECT RECORDS AND BUILD A PROFILE", 20 * "-")
    print ("1. Set reference sequence for each feature")
    print ("2. Select records to build profile for")
    print ("3. Select features of selected records")
    print ("4. Print summary of gene records")
    print ("5. Build alignment and profile")
    print ("6. Use profile to search genome")
    print ("7. Return to main menu")
    print (77 * "-")
    choice = input(" >>  ")
    choice = profile_dict[choice]
    exec_menu(choice)
    return

def set_reference():
    global referenceSeqs
    referenceSeqs = {}
    feature = input("Which feature do you want to set a reference sequence for? Available features are " + str([i for i in user_annotations]))
    ref_sequence = input ("Which record should be used as the reference sequence? Available records are" + str([i for i in seqDict if i in feature_to_record[feature]]))
    referenceSeqs[feature] = ref_sequence
    print (ref_sequence, "has been set as the reference sequence for", feature)
    profile_menu()

def select_records():
    print ("\n")
    print (30 * "-", "SELECT RECORDS", 30 * "-")
    print ("1. Only allow unique strains")
    print ("2. Only allow unique species")
    print ("3. Return to profile menu")
    print (67 * "-")
    choice = input(" >>  ")
    choice = record_dict[choice]
    exec_menu(choice)
    return

    profile_menu()


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

    profile_menu()


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

    profile_menu()

def getAlignmentScore(seq1, seq2):
    print (seq1)
    print (seq2)

    alignment_score = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5, score_only=True)
    return alignment_score



def select_features():
    global selected_records, selected_features
    required_features = input("Only select this feature ")
    if required_features not in user_annotations:
        print ("\n *** ERROR: no records have this feature. \n")
    else:
        selected_features.clear()
        selected_features.add(required_features)
    profile_menu()

def feature_summary():
    global selected_records, selected_features
    if len(selected_features) == 0:
        selected_features = user_annotations

    if len(selected_records) == 0:
        selected_records = seqDict
    print ("Record summary")
    print (selected_records.keys())
    print ("Feature summary")
    print (selected_features)
    profile_menu()

def build_profile():
    global selected_records, selected_features
    align_list = []


    if len(selected_features) == 0:
        selected_features = user_annotations


    if len(selected_records) == 0:
        selected_records = seqDict

    for feature in selected_features:
        # print (feature)
        for record in selected_records:
            if feature in selected_records[record].annotations.keys():

                align_record = SeqRecord(Seq(selected_records[record].annotations[feature], generic_protein), id=record + "_" + feature)
                align_list.append(align_record)

        # Write sequences to FASTA

        SeqIO.write(align_list, "align.fasta", "fasta")



        # Align using ClustalOmega
        muscle_cline = MuscleCommandline(input="align.fasta")
        child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE ,stderr=subprocess.PIPE, universal_newlines=True, shell = (sys.platform!="win32"))
        alignment = AlignIO.read(child.stdout, "fasta")
        AlignIO.write(alignment, "align.aln", "fasta")
        subprocess.call(["hmmbuild", "profile.hmm", "align.aln"])
        # print (alignment)



    profile_menu()




def summary_menu():
    print (30 * "-", " PRINT SUMMARY", 30 * "-")
    print ("1. Print summary of all records")
    print ("2. Print summary of only selected records")
    print ("3. Return to main menu")
    print (67 * "-")
    choice = input(" >>  ")
    choice = summary_dict[choice]
    exec_menu(choice)
    return




def summary_all():
    print ("Here are all the loaded records: \n ")
    printRecords(seqDict, user_annotations)

    if len(unmappable) > 0:
        print ("The following IDs couldn't be mapped to a GenBank record: ")
        for record in unmappable:
            print (record)
        print ("\n")
    summary_menu()


def summary_selected():
    global selected_featur3es

    if len(selected_records) == 0:
        print ("You haven't selected any records")
        summary_menu()
        return

    else:
        print ("Here are all the selected records: ")
        printRecords(selected_records, selected_features)

        if len(unmappable) > 0:
            print ("The following IDs couldn't be mapped to a GenBank record: ")
            for record in unmappable:
                print (record)
            print ("\n")
    summary_menu()


def printRecords(seqRecords, selected_features):
    if len(selected_features) == 0:
        selected_features = user_annotations



    for record in seqRecords:
        feature_present = any(i in selected_features for i in seqRecords[record].annotations.keys())

        if feature_present:
            print (seqRecords[record].id, "-", seqRecords[record].annotations['organism'])
            for key in seqRecords[record].annotations.keys():
                if key in selected_features:
                    print ("Has feature", key)
                    if 'Overlap' in key:
                        if key:
                            print (key.split()[0], "features are fused")
                        else:
                            print (key.split()[0], "features are split")
                            print ("And the distance between them is ", seqRecords[record].annotations[key.split()[0] + "Distance"])
        print ("\n")

def search_genome():
    child = subprocess.Popen("", stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                             shell=(sys.platform != "win32"))
    subprocess.call(["hmmsearch", "profile.hmm", "genome.fa"])

    profile_menu()


# Back to main menu
def back():
    menu_actions['main_menu']()


# Exit program
def exit():
    sys.exit()


# =======================
#    MENUS DEFINITIONS
# =======================
# Menu definition
menu_actions = {
    'main_menu': main_menu,
    'back': back,
    '1': map_to_record,
    '2': annotate_menu,
    '3': profile_menu,
    '4': summary_menu,
    '5': exit,
    '6': overlap_menu,
    '10' : select_records,
    '11' : select_features,
    '12' : feature_summary,
    '13' : build_profile,
    '14' : unique_strains,
    '15' : unique_species,
    '16' : summary_all,
    '17' : summary_selected,
    '18' : set_reference,
    '19' : search_genome

}

# =======================
#      MAIN PROGRAM
# =======================

# Main Program
if __name__ == "__main__":


    # Load default values
    # defaultValues("files/YenA1_tiny.fasta", "A1", "files/YenA2_tiny.fasta", "A2", )
    # print ("\n *** Default values loaded *** \n")
    # Launch main menu
    main_menu()