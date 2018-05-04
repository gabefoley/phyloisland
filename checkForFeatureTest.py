import checkForFeature
from Bio.SearchIO import HmmerIO
from Bio import SearchIO

# checkForFeature.read_hmmer_results("tmp/FN667742.1GNTPI_backward_0_hmmsearch_results.fasta")

# hmm_parser = HmmerIO.Hmmer3TabParser("tmp/FN667742.1GNTPI_backward_0_translated_genome.fasta")

hmmer_results = SearchIO.parse("tmp/FN667742.1GNTPI_backward_0_translated_genome.fasta", "hmmer3-text")


print (hmmer_results)
# for x in hmm_parser:
#     print (x)