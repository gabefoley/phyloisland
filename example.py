from Bio import SeqIO

# Load in a file with a single sequence
single = SeqIO.parse("files/YenA1_single.fasta", "fasta")

# Load in a file with multiple sequences
records = SeqIO.parse("files/YenA1_test.fasta", "fasta")


# Function which takes a letter and a sequence and counts how many times that letter occurs in the sequence
def count_chars(letter, seq):
    count = seq.count(letter)

    # This part just works out if we need to add a plural to the letter
    plural = ""
    if count != 1:
        plural = "'s"

    # Return a string saying what we found
    return("We found %s %s%s in here" % (count, letter, plural))






# Print out the full record from the single sequence file
for record in single:
    print (record)
    print()

# Print out just the record name and first 10 characters for the file with multiple sequences
for record in records:
    print (record.name, record.seq[0:20])
    print (count_chars("Q", record.seq[0:20]))
    print ()



