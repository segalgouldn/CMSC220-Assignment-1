# Noah Segal-Gould <ns2349@bard.edu>
# 2018-02-01
# CMSC 220
# Lab Assignment 1
# My code is based in part on work done collectively in class.

# I added some more sequences for testing purposes.
sequence_1 = "CTTTCGCAATGATAC"
sequence_2 = "GAAAGCGTTACTATG"
sequence_3 = "UAGGCCAUCCGUAAC"
sequence_4 = "CTTTCGCAAZZZZZZ"


# I modified the validation step.
def ligase(seq_1, seq_2):
    sum_sequence = seq_1 + seq_2
    if all(n in "ATCG" for n in sum_sequence):
        return sum_sequence
    else:
        raise Exception("DNA sequence invalid.")
    

# The function accepts both RNA and DNA, and also replaces non-nucleotides with "N"
def reverseTransDNA(seq):
    seq = seq.upper()
    if "U" not in seq:
        rna_sequence = seq.replace("T", "U")
        translation_dictionary = {"A": "U", "U": "A", "C": "G", "G": "C"}
        return seq, "".join([translation_dictionary[n] if n in "AUCG" else "N" for n in reversed(rna_sequence)])
    else:
        dna_sequence = seq.replace("U", "T")
        translation_dictionary = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return seq, "".join([translation_dictionary[n] if n in "ATCG" else "N" for n in reversed(dna_sequence)])

    
# This function splits sequences in two parts at the chosen index.
def nuclease(seq, cutpos):
    cut_1 = seq[:cutpos]
    cut_2 = seq[cutpos:]
    return cut_1, cut_2


# Tests
print("*"*10 + " TESTS " + "*"*10)
print(sequence_1 + " + " + sequence_2 + " => " + ligase(sequence_1, sequence_2))
print("DNA to RNA of " + sequence_1 + ": " + str(reverseTransDNA(sequence_1)))
print("DNA to RNA of " + sequence_2 + ": " + str(reverseTransDNA(sequence_2)))
print("RNA to DNA of " + sequence_3 + ": " + str(reverseTransDNA(sequence_3)))
print("DNA to RNA of " + sequence_4 + " (with errors): " + str(reverseTransDNA(sequence_4)))
print("Nuclease split on " + sequence_1 + ": " + str(nuclease(sequence_1, 3)))
print("Nuclease split on " + sequence_2 + ": " + str(nuclease(sequence_2, 3)))
