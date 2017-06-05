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
from Bio.pairwise2 import format_alignment



# from BioSQL import BioSeqDatabase
# server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi",
#                      passwd = "", host = "localhost", db="phyloisland")
# db = server.new_database("phylomain", description="Just for testing")
# server.commit() #On Biopython 1.49 or older, server.adaptor.commit()

from Bio import pairwise2
alignments = pairwise2.align.globalms("AAAAGGGGCGA", "AAGGT", 2, -1, -.5, -.1)
print(format_alignment(*alignments[0]))