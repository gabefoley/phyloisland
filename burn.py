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



seq1 = Seq("AAPPNAGAGAGA", generic_protein)
seq2 = Seq("AAPVVVVWWWWPNAGAGANNNNGA", generic_protein)


alignments = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5, score_only=True)
print (alignments)
print (type(alignments))
