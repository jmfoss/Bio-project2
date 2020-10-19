def main():
    n_gene = get_gene(28274, 29533)


# param: start and end index of dna segment
# return: string of rna segment 3'-5'
def get_gene(begin, end):
    sequence = open("sequence.txt", "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]


# Causes the program to run main when script is ran
if __name__ == '__main__':
    main()