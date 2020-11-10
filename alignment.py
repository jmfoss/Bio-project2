from Bio import pairwise2
from Bio.pairwise2 import format_alignment


class Alignment:
    seqA = str()
    seqB = str()
    score = int()
    start = int()
    end = int()


# Dictionary of codons and amino acids
table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


# Given two sequences, returns alignment object
def get_global_alignment(seqA, seqB):
    # Number of columns and rows
    col = ' ' + seqA
    row = ' ' + seqB

    # 0: align
    # 1 : Gap top
    # 2 : Gap low

    # Creates matrix for finding alignment
    matrix = [[(0, 0) for i in range(len(row))] for j in range(len(col))]

    for indexB in range(len(row)):
        for indexA in range(len(col)):
            # Checks if both are empty
            if col[indexA] == ' ' and row[indexB] == ' ':
                break
            # Checks if seqA is empty
            elif col[indexA] == ' ':
                matrix[indexA][indexB] = (1, matrix[indexA][indexB - 1][1] - 2)
            # Checks if seqB iis empty
            elif row[indexB] == ' ':
                matrix[indexA][indexB] = (2, matrix[indexA - 1][indexB][1] - 2)
            else:
                # Check align
                if col[indexA] == row[indexB]:
                    # if match
                    matrix[indexA][indexB] = (0, matrix[indexA - 1][indexB - 1][1] + 1)
                else:
                    # if mismatch
                    matrix[indexA][indexB] = (0, matrix[indexA - 1][indexB - 1][1] - 1)

                # Check top gap
                local_score = matrix[indexA][indexB - 1][1] - 2
                matrix[indexA][indexB] = (1, local_score) if local_score > matrix[indexA][indexB][1] else \
                    matrix[indexA][indexB]

                # Check low gap
                local_score = matrix[indexA - 1][indexB][1] - 2
                matrix[indexA][indexB] = (2, local_score) if local_score > matrix[indexA][indexB][1] else \
                    matrix[indexA][indexB]
    # Used for iterating through matrix backwards
    indexA = len(col) - 1
    indexB = len(row) - 1

    # Creates alignment object
    alignment = Alignment()
    # Assigns score
    alignment.score = matrix[indexA][indexB][1]

    # Finds alignment backwards
    while indexA != 0 and indexB != 0:
        if matrix[indexA][indexB][0] == 0:
            alignment.seqA += col[indexA]
            alignment.seqB += row[indexB]
            indexA -= 1
            indexB -= 1
        elif matrix[indexA][indexB][0] == 1:
            alignment.seqA += '-'
            alignment.seqB += row[indexB]
            indexB -= 1
        else:
            alignment.seqA += col[indexA]
            alignment.seqB += '-'
            indexA -= 1

    # Reverses sequences
    alignment.seqA = alignment.seqA[::-1]
    alignment.seqB = alignment.seqB[::-1]

    # Assigns start and end points
    alignment.start = 0
    alignment.end = len(alignment.seqA)
    return alignment


# Given an alignment object
# Return number of insertion/deletions
def countIndels(alignment):
    insertions = 0
    deletions = 0
    prev = False
    for nuc in alignment.seqA:      
        if nuc == '-':
            if prev == False:
                # A gap means there was a deletion in one sequence, but an insertion
                # for the other, so increment both.
                insertions += 1
                prev = True
            else:
                prev = False

    prev = False
    for nuc in alignment.seqB:            
        if nuc == '-':
            if prev == False:
                # A gap means there was a deletion in one sequence, but an insertion
                # for the other, so increment both.
                deletions += 1
                prev = True
            else:
                prev = False

    return insertions + deletions


def main():
    # Loading genes
    n_mers_gene = get_gene(28566, 29807, "MERS_sequence.txt")
    n_covid_gene = get_gene(28274, 29533, "COVID_sequence.txt")

    # Aligning genes
    global_alignment = get_global_alignment(n_mers_gene, n_covid_gene)
    middle = str()
    for a, b in zip(global_alignment.seqA, global_alignment.seqB):
        if a == b:
            middle += '|'
        elif a == '-' or b == '-':
            middle += ' '
        else:
            middle += '.'
    print("Our alignment")
    print(global_alignment.seqA)
    print(middle)
    print(global_alignment.seqB)

    mutations = find_mutations(global_alignment)
    indel = countIndels(global_alignment)

    print("  Score=" + str(global_alignment.score))
    print("\nindel: " + str(indel))
    print("non-syn mutations: " + str(mutations[0]))
    print("syn mutations: " + str(mutations[1]))

    print("\nBio Python Alignment")
    global_alignments = pairwise2.align.globalms(n_mers_gene, n_covid_gene, 1, -1, -2, -2)
    print(format_alignment(*global_alignments[0]))
    mutations = find_mutations(global_alignments[0])
    indel = countIndels(global_alignments[0])
    print("indel: " + str(indel))
    print("non-syn mutations: " + str(mutations[0]))
    print("syn mutations: " + str(mutations[1]))


# Given alignment object
# Returns tuple of mutations (non_syn, syn)
def find_mutations(alignment):
    # First: non_syn
    # Second: syn
    mutations = [0, 0]
    for codon in range(alignment.start, alignment.end, 3):
        codonA = alignment.seqA[codon: codon + 3]
        codonB = alignment.seqB[codon: codon + 3]
        if codonA != codonB and codonA in table and codonB in table:
            if table[codonA] != table[codonB]:
                mutations[0] += 1
            else:
                mutations[1] += 1
    return mutations


# param: start and end index of dna segment
# return: string of rna segment 3'-5'
def get_gene(begin, end, filename):
    sequence = open(filename, "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]


# Causes the program to run main when script is ran
if __name__ == '__main__':
    main()
