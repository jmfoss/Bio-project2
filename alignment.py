from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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


def get_global_alignment(seqA, seqB):
    col = ' ' + seqA
    row = ' ' + seqB
    # 0: align
    # 1 : Gap top
    # 2 : Gap low
    matrix = [[(0, 0) for i in range(len(row))] for j in range(len(col))]
    for indexB in range(len(row)):
        for indexA in range(len(col)):
            if col[indexA] == ' ' and row[indexB] == ' ':
                break
            elif col[indexA] == ' ':
                matrix[indexA][indexB] = (1, matrix[indexA][indexB - 1][1] - 2)
            elif row[indexB] == ' ':
                matrix[indexA][indexB] = (2, matrix[indexA - 1][indexB][1] - 2)
            else:
                # Check align
                if col[indexA] == row[indexB]:
                    matrix[indexA][indexB] = (0, matrix[indexA - 1][indexB - 1][1] + 1)
                else:
                    matrix[indexA][indexB] = (0, matrix[indexA - 1][indexB - 1][1] - 1)

                # Check top gap
                local_score = matrix[indexA][indexB - 1][1] - 2
                matrix[indexA][indexB] = (1, local_score) if local_score > matrix[indexA][indexB][1] else matrix[indexA][indexB]

                # Check low gap
                local_score = matrix[indexA - 1][indexB][1] - 2
                matrix[indexA][indexB] = (2, local_score) if local_score > matrix[indexA][indexB][1] else matrix[indexA][indexB]
    indexA = len(col) - 1
    indexB = len(row) - 1
    alignment = dict()
    alignment['seqA'] = str()
    alignment['seqB'] = str()
    print(matrix)
    alignment['score'] = matrix[indexA][indexB][1]
    while indexA != 0 and indexB != 0:
        if matrix[indexA][indexB][0] == 0:
            alignment['seqA'] += col[indexA]
            alignment['seqB'] += row[indexB]
            indexA -= 1
            indexB -= 1
        elif matrix[indexA][indexB][0] == 1:
            alignment['seqA'] += '-'
            alignment['seqB'] += row[indexB]
            indexB -= 1
        else:
            alignment['seqA'] += col[indexA]
            alignment['seqB'] += '-'
            indexA -= 1

    alignment['seqA'] = alignment['seqA'][::-1]
    alignment['seqB'] = alignment['seqB'][::-1]
    alignment['start'] = 0
    alignment['end'] = len(alignment['seqA'])
    return alignment

def countIndels(seqA, seqB):
    insertions = 0
    deletions = 0
    for nuc in seqA:
        if nuc == '-':
            # A gap means there was a deletion in one sequence, but an insertion
            # for the other, so increment both.
            insertions += 1
            deletions += 1
            
    for nuc in seqB:
        if nuc == '-':
            # A gap means there was a deletion in one sequence, but an insertion
            # for the other, so increment both.
            insertions += 1
            deletions += 1
    
    return insertions + deletions

def main():
    # Loading genes
    n_mers_gene = get_gene(28566, 29807, "MERS_sequence.txt")
    n_covid_gene = get_gene(28274, 29533, "COVID_sequence.txt")

    # Aligning genes
    global_alignment = get_global_alignment(n_mers_gene, n_covid_gene)

    # Mutation count
    non_syn = 0
    syn = 0
    indel = 0
    # Printing alignment
    for codon in range(global_alignment['start'], global_alignment['end'], 3):
        codonA = global_alignment['seqA'][codon: codon + 3]
        codonB = global_alignment['seqB'][codon: codon + 3]
        if codonA != codonB and codonA in table and codonB in table:
            if table[codonA] != table[codonB]:
                non_syn += 1
            else:
                syn += 1
        indel = countIndels(global_alignment['seqA'], global_alignment['seqB'])

    print(global_alignment)
    print(indel)


# param: start and end index of dna segment
# return: string of rna segment 3'-5'
def get_gene(begin, end, filename):
    sequence = open(filename, "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]


# Causes the program to run main when script is ran
if __name__ == '__main__':
    main()
