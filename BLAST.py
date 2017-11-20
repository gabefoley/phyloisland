# First we must import the Seq module
from Bio import SeqIO

# Load in a file with a single sequence
my_seq = SeqIO.parse("files/YenA1_single.fasta", "fasta")

# Now we must specify the correct alphabet for the sequence
from Bio.Alphabet import IUPAC

# pre running checks
# print your fasta to verify it is the correct one
for seq in my_seq:
    print(seq.id)
    print(repr(seq.seq))
    print(len(seq))

# Time to use local BLAST
from Bio.Blast.Applications import NcbiblastxCommandline
blastn_cline = NcbiblastxCommandline(query="my_seq", db="nr", evalue=0.001,
                                      outfmt=5, out="my_seq_blast_results.xml")
blastx_cline = NcbiblastxCommandline(cmd='blastn', out='opuntia.xml', outfmt=5, query='opuntia.fasta',
db='nr', evalue=0.001)
print(blastn_cline)
# blastn -out my_seq_blast_results.xml -outfmt 5 -query opuntia.fasta -db nr -evalue 0.001
# stdout, stderr = blastx_cline()

# create a handle to blast and open the file for viewing
result_handle = open("my_seq_blast_results.xml")

# Now you need to parse the blast output
from Bio.Blast import NCBIXML
blast_records = NCBIXML.read(result_handle)

# This will help to step through large BLAST outputs
for blast_record in blast_records:
     # Do something with blast_record
    blast_records = list(blast_records)

# Check out some of you BLAST result info
E_VALUE_THRESH = 0.04 #you can adjust this
for alignment in blast_record.alignments:
     for hsp in alignment.hsps:
         if hsp.expect < E_VALUE_THRESH:
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')
