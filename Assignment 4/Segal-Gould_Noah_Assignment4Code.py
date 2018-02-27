
# Noah Segal-Gould <ns2349@bard.edu>
# 2018-02-27
# CMSC 220
# Origin of Replication Exercise
# My code is based in part on work done collectively in class as well as code provided on Wikipedia and the course Moodle 2.

# Import libraries for managing regular expressions and counting k-mers
import re
from collections import Counter

# Create a dictionary to keep track of k-mers of each size (3, 6, and 9) in each reading frame (1-6)
def find_reading_frames_k_mers(sequence):
    sequence = sequence.upper().replace("\n", "") # Capitalize nucleotides and remove newlines
    translation_table = str.maketrans("ATCG", "TAGC")
    return {1: {3: Counter(re.findall(".{3}", sequence)).most_common(4), 
                6: Counter(re.findall(".{6}", sequence)).most_common(4), 
                9: Counter(re.findall(".{9}", sequence)).most_common(4)},
            2: {3: Counter(re.findall(".{3}", sequence[1:-2])).most_common(4), 
                6: Counter(re.findall(".{6}", sequence[1:-2])).most_common(4), 
                9: Counter(re.findall(".{9}", sequence[1:-2])).most_common(4)},
            3: {3: Counter(re.findall(".{3}", sequence[2:-1])).most_common(4), 
                6: Counter(re.findall(".{6}", sequence[2:-1])).most_common(4), 
                9: Counter(re.findall(".{9}", sequence[2:-1])).most_common(4)},
            4: {3: Counter(re.findall(".{3}", sequence.translate(translation_table))[::-1]).most_common(4), 
                6: Counter(re.findall(".{6}", sequence.translate(translation_table))[::-1]).most_common(4), 
                9: Counter(re.findall(".{9}", sequence.translate(translation_table))[::-1]).most_common(4)},
            5: {3: Counter(re.findall(".{3}", sequence[2:-1].translate(translation_table))[::-1]).most_common(4), 
                6: Counter(re.findall(".{6}", sequence[2:-1].translate(translation_table))[::-1]).most_common(4), 
                9: Counter(re.findall(".{9}", sequence[2:-1].translate(translation_table))[::-1]).most_common(4)},
            6: {3: Counter(re.findall(".{3}", sequence[1:-2].translate(translation_table))[::-1]).most_common(4), 
                6: Counter(re.findall(".{6}", sequence[1:-2].translate(translation_table))[::-1]).most_common(4), 
                9: Counter(re.findall(".{9}", sequence[1:-2].translate(translation_table))[::-1]).most_common(4)}}

# Returns a list of k-mers based on some sequence and a value "k"
# Adapted from Wikipedia
def find_k_mers(sequence, k):
    k_mers = []
    n = len(sequence)
    for i in range(0, n-k+1):
        k_mers.append(sequence[i:i+k])
    return k_mers

# Ignores the other 5 reading frames and puts the 3, 6, and 9-mers of just the first reading frame in a dictionary
def find_non_reading_frames_k_mers(sequence):
    sequence = sequence.upper().replace("\n", "") # Capitalize nucleotides and remove newlines
    results = {1: {}}
    for j in [3, 6, 9]:
        results[1][j] = Counter(find_k_mers(sequence, j)).most_common(4)
    return results

# Print the dictionaries as tables
def print_table(sequence, use_reading_frames=True):
    if use_reading_frames:
        for reading_frame, k_mers_dict in find_reading_frames_k_mers(sequence).items():
            print("{0}|3|{1}         {2}         {3}         {4}\n |6|{5}      {6}      {7}      {8}\n |9|{9}   {10}   {11}   {12}"
                  .format(reading_frame, 
                          k_mers_dict[3][0][0], k_mers_dict[3][1][0], k_mers_dict[3][2][0], k_mers_dict[3][3][0],
                          k_mers_dict[6][0][0], k_mers_dict[6][1][0], k_mers_dict[6][2][0], k_mers_dict[6][3][0],
                          k_mers_dict[9][0][0], k_mers_dict[9][1][0], k_mers_dict[9][2][0], k_mers_dict[9][3][0]))
            print(" " + "="*48)
    else:
        for reading_frame, k_mers_dict in find_non_reading_frames_k_mers(sequence).items():
            print("{0}|3|{1}         {2}         {3}         {4}\n |6|{5}      {6}      {7}      {8}\n |9|{9}   {10}   {11}   {12}"
                  .format(reading_frame, 
                          k_mers_dict[3][0][0], k_mers_dict[3][1][0], k_mers_dict[3][2][0], k_mers_dict[3][3][0],
                          k_mers_dict[6][0][0], k_mers_dict[6][1][0], k_mers_dict[6][2][0], k_mers_dict[6][3][0],
                          k_mers_dict[9][0][0], k_mers_dict[9][1][0], k_mers_dict[9][2][0], k_mers_dict[9][3][0]))
            print(" " + "="*48)

# Main function to run tests
def main():
    sequence_vibrio_cholerae = """atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaac
    ctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttct
    tggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggatt
    acgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaa
    gccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaag
    atcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaat
    gatcaagctgctgctcttgatcatcgtttc"""
    
    sequence_thermotoga_petrophila = """aactctatacctcctttttgtcgaatttgtgtgatttatagagaaaatcttatt
    aactgaaactaaaatggtaggtttggtggtaggttttgtgtacattttgtagtatctgatttttaattacataccgtatattgtattaaa
    ttgacgaacaattgcatggaattgaatatatgcaaaacaaacctaccaccaaactctgtattgaccattttaggacaacttcagggtggt
    aggtttctgaagctctcatcaatagactattttagtctttacaaacaatattaccgttcagattcaagattctacaacgctgttttaatg
    ggcgttgcagaaaacttaccacctaaaatccagtatccaagccgatttcagagaaacctaccacttacctaccacttacctaccacccgg
    gtggtaagttgcagacattattaaaaacctcatcagaagcttgttcaaaaatttcaatactcgaaacctaccacctgcgtcccctattat
    ttactactactaataatagcagtataattgatctgaaaagaggtggtaaaaaa"""
    
    vibrio_cholerae_header_1 = "K-Mers for Vibrio Cholerae, with reading frames:"
    print(vibrio_cholerae_header_1)
    print("-"*len(vibrio_cholerae_header_1))
    print_table(sequence_vibrio_cholerae, use_reading_frames=True)
    
    print()
    
    vibrio_cholerae_header_2 = "K-Mers for Vibrio Cholerae, without reading frames:"
    print(vibrio_cholerae_header_2)
    print("-"*len(vibrio_cholerae_header_2))
    print_table(sequence_vibrio_cholerae, use_reading_frames=False)
    
    print()
    
    thermotoga_petrophila_header_1 = "K-Mers for Thermotoga Petrophila, with reading frames:"
    print(thermotoga_petrophila_header_1)
    print("-"*len(thermotoga_petrophila_header_1))
    print_table(sequence_thermotoga_petrophila, use_reading_frames=True)
    
    print()
    
    thermotoga_petrophila_header_2 = "K-Mers for Thermotoga Petrophila, without reading frames:"
    print(thermotoga_petrophila_header_2)
    print("-"*len(thermotoga_petrophila_header_2))
    print_table(sequence_thermotoga_petrophila, use_reading_frames=False)

# Run the main function
main()