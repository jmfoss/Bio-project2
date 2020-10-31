from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def main():
    n_mers_gene = get_gene(28566, 29807, "MERS_sequence.txt")
    n_covid_gene = get_gene(28274, 29533, "COVID_sequence.txt")
    global_alignments = pairwise2.align.globalms(n_mers_gene, n_covid_gene, 1, -1, -2, -2)
    local_alignments = pairwise2.align.localms(n_mers_gene, n_covid_gene, 1, -1, -2, -2)


# param: start and end index of dna segment
# return: string of rna segment 3'-5'
def get_gene(begin, end, filename):
    sequence = open(filename, "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]


# Causes the program to run main when script is ran
if __name__ == '__main__':
    main()
