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
            print(record)
        print("\n")
    summary_menu()


def summary_selected():
    global selected_featur3es

    if len(selected_records) == 0:
        print("You haven't selected any records")
        summary_menu()
        return

    else:
        print("Here are all the selected records: ")
        printRecords(selected_records, selected_features)

        if len(unmappable) > 0:
            print("The following IDs couldn't be mapped to a GenBank record: ")
            for record in unmappable:
                print(record)
            print("\n")
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