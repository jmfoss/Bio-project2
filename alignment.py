def main():
    n_mers_gene = get_gene(28274, 29533, "MERS_sequence.txt")
    n_covid_gene = get_gene(28274, 29533, "COVID_sequence.txt")

# param: start and end index of dna segment
# return: string of rna segment 3'-5'
def get_gene(begin, end, filename):
    sequence = open(filename, "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]


# Causes the program to run main when script is ran
if __name__ == '__main__':
    main()