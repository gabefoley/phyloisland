from Bio import Entrez, SeqIO, GenBank, AlignIO, pairwise2
Entrez.email = "gabriel.foley@uqconnect.edu.au"



handle = Entrez.efetch(db='nucleotide', rettype='gb', id='CP010029.1', retmode='text')
protein_records = SeqIO.parse(handle, "gb")

for record in protein_records:
    print(record.description)
    for feature in record.features:
        if 'translation' in feature.qualifiers:
            if (feature.qualifiers['translation'][0].startswith("MSNSIEAKLQEDLRDALVDYYLGQIVP")):
                print (feature)
                print (feature.qualifiers['translation'][0])
print('*******************')