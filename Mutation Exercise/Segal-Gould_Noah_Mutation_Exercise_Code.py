
# Noah Segal-Gould <ns2349@bard.edu>
# 2018-02-27
# CMSC 220
# Mutation Exercise
# My code is based in part on work done collectively in class as well as code provided on the course Moodle 2.

# Import libraries for managing random indices, random list elements, and regular expressions
from random import choice, randint
import re

# A list of possible nucleotides from which random selections will be taken
nucleotides = list("ATCG")

# Dictionary taken from Moodle 2
standard_genetic_code = {'UUU':'Phe', 'UUC':'Phe', 'UCU':'Ser', 'UCC':'Ser',
                         'UAU':'Tyr', 'UAC':'Tyr', 'UGU':'Cys', 'UGC':'Cys',
                         'UUA':'Leu', 'UCA':'Ser', 'UAA':'Stop','UGA':'Stop',
                         'UUG':'Leu', 'UCG':'Ser', 'UAG':'Stop','UGG':'Trp',
                         'CUU':'Leu', 'CUC':'Leu', 'CCU':'Pro', 'CCC':'Pro',
                         'CAU':'His', 'CAC':'His', 'CGU':'Arg', 'CGC':'Arg',
                         'CUA':'Leu', 'CUG':'Leu', 'CCA':'Pro', 'CCG':'Pro',
                         'CAA':'Gln', 'CAG':'Gln', 'CGA':'Arg', 'CGG':'Arg',
                         'AUU':'Ile', 'AUC':'Ile', 'ACU':'Thr', 'ACC':'Thr',
                         'AAU':'Asn', 'AAC':'Asn', 'AGU':'Ser', 'AGC':'Ser',
                         'AUA':'Ile', 'ACA':'Thr', 'AAA':'Lys', 'AGA':'Arg',
                         'AUG':'Met', 'ACG':'Thr', 'AAG':'Lys', 'AGG':'Arg',
                         'GUU':'Val', 'GUC':'Val', 'GCU':'Ala', 'GCC':'Ala',
                         'GAU':'Asp', 'GAC':'Asp', 'GGU':'Gly', 'GGC':'Gly',
                         'GUA':'Val', 'GUG':'Val', 'GCA':'Ala', 'GCG':'Ala', 
                         'GAA':'Glu', 'GAG':'Glu', 'GGA':'Gly', 'GGG':'Gly'}

# Dictionary taken from Moodle 2
amino_acid_code = {'Ala':'A', 'Cys':'C', 'Asp':'D', 'Glu':'E',
                   'Phe':'F', 'Gly':'G', 'His':'H', 'Ile':'I',
                   'Lys':'K', 'Leu':'L', 'Met':'M','Asn':'N',
                   'Pro':'P', 'Gln':'Q', 'Arg':'R','Ser':'S',
                   'Thr':'T', 'Val':'V', 'Trp':'W', 'Tyr':'Y',
                   'Stop':'_'}

# Randomly generate nucleotide sequences which are ten characters long
def generate():
    return "".join([choice(nucleotides) for n in range(10)])

# Based on user choice, mutate some sequence in one of three ways
def mutate(sequence, mutation_type):
    index = randint(0, len(sequence)-1)
    if mutation_type == "point":
        mutation = sequence[:index] + choice(nucleotides) + sequence[index + 1:]
        if mutation == sequence:
            mutation = mutate(sequence, "point") # Recursive call; should keep trying until it's successfully mutated
    elif mutation_type == "deletion":
        mutation = sequence[:index] + sequence[index + 1:]
    elif mutation_type == "insertion":
        mutation = sequence[:index] + choice(nucleotides) + sequence[index:]
    else: 
        raise Exception("Invalid Mutation Type: {}".format(mutation_type))
    return mutation

# Identify which kind of mutation occurred
def detect(sequence, mutation):
    # If they are the same, there was no mutation
    if sequence != mutation:
        # Find the amino acid strings of both sequences
        mutated_amino_acids = "".join([amino_acid_code[standard_genetic_code[codon]] 
                                       for codon in re.findall(".{3}", mutation.replace("T", "U"))])
        sequence_amino_acids = "".join([amino_acid_code[standard_genetic_code[codon]] 
                                       for codon in re.findall(".{3}", sequence.replace("T", "U"))])
        if mutated_amino_acids == sequence_amino_acids:
            return "Going from {0} to {1} there was a silent mutation".format(sequence, mutation)
        elif len(mutation) > len(sequence) and mutated_amino_acids != sequence_amino_acids:
            return "Going from {0} to {1} there was a frameshift insertion mutation".format(sequence, mutation)
        elif len(mutation) < len(sequence) and mutated_amino_acids != sequence_amino_acids:
            return "Going from {0} to {1} there was a frameshift deletion mutation".format(sequence, mutation)
        elif mutated_amino_acids != sequence_amino_acids:
            return "Going from {0} to {1} there was a missense mutation".format(sequence, mutation)
        elif "_" in mutation and "_" not in sequence:
            return "Going from {0} to {1} there was a nonsense mutation".format(sequence, mutation)
        else:
            raise Exception("No mutation occurred going from {0} to {1}".format(sequence, mutation))
    else:
        raise Exception("No mutation occurred going from {0} to {1}".format(sequence, mutation))

# Main function for running tests
def main():
    sequence = "ATCGATCGAT"
    generated = generate()
    point_mutated_sequence = mutate(sequence, "point")
    point_mutated_generated_sequence = mutate(generated, "point")
    deletion_mutated_sequence = mutate(sequence, "deletion")
    deletion_mutated_generated_sequence = mutate(generated, "deletion")
    insertion_mutated_sequence = mutate(sequence, "insertion")
    insertion_mutated_generated_sequence = mutate(generated, "insertion")
    print("="*20 + " Testing Mutations " + "="*20)
    print("The sequences are {0} & {1}".format(sequence, generated))
    print("The sequence \"{0}\" undergoes a point mutation and becomes \"{1}\"".format(sequence, point_mutated_sequence))
    print("The sequence \"{0}\" undergoes a deletion mutation and becomes \"{1}\"".format(sequence, deletion_mutated_sequence))
    print("The sequence \"{0}\" undergoes an insertion mutation and becomes \"{1}\"".format(sequence, insertion_mutated_sequence))
    print("The sequence \"{0}\" undergoes a point mutation and becomes \"{1}\"".format(generated, point_mutated_generated_sequence))
    print("The sequence \"{0}\" undergoes a deletion mutation and becomes \"{1}\"".format(generated, deletion_mutated_generated_sequence))
    print("The sequence \"{0}\" undergoes an insertion mutation and becomes \"{1}\"".format(generated, insertion_mutated_generated_sequence))
    print("="*20 + " Testing Mutation Detection " + "="*20)
    print(detect(sequence, point_mutated_sequence))
    print(detect(sequence, deletion_mutated_sequence))
    print(detect(sequence, insertion_mutated_sequence))
    print(detect(generated, point_mutated_generated_sequence))
    print(detect(generated, deletion_mutated_generated_sequence))
    print(detect(generated, insertion_mutated_generated_sequence))

# Run the main function
main()