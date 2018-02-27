# Noah Segal-Gould <ns2349@bard.edu>
# 2018-02-17
# CMSC 220
# Lab Assignment 3
# My code is based in part on work done collectively in class as well as code provided on the course Moodle 2.

# Import regular expressions library
import re

class Sequence: # Create a Sequence class
    def __init__(self, filename): # Create constructor
        self.filename = filename # Set filename as attribute
        with open(filename) as f: # Open the file using the filename
            # Dictionaries provided on Moodle for converting from nucleotides to amino acids
            self.standard_genetic_code = {'UUU':'Phe', 'UUC':'Phe', 'UCU':'Ser', 'UCC':'Ser',
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
            self.amino_acid_code = {'Ala':'A', 'Cys':'C', 'Asp':'D', 'Glu':'E',
                                    'Phe':'F', 'Gly':'G', 'His':'H', 'Ile':'I',
                                    'Lys':'K', 'Leu':'L', 'Met':'M','Asn':'N',
                                    'Pro':'P', 'Gln':'Q', 'Arg':'R','Ser':'S',
                                    'Thr':'T', 'Val':'V', 'Trp':'W', 'Tyr':'Y',
                                    'Stop':'_'}
            # Remove unnecessary data from file and clean up nucleotide text
            self.DNA = "".join([line for line in f.read().upper().split("\n") if ">" not in line])
            # Prepare for transcription 
            self.RNA = self.DNA.replace("T", "U")
            # Calculate some useful statistics
            self.nucleotide_totals = {"A": self.DNA.count("A"), 
                                      "T": self.DNA.count("T"), 
                                      "C": self.DNA.count("C"), 
                                      "G": self.DNA.count("G"), 
                                      "N": len(re.findall("[^ATCGatcg]", self.DNA))}
            # Only make reading frames if there aren't any non-nucleotides in the sequence
            if self.nucleotide_totals["N"] == 0:
                self.translation_table = str.maketrans("AUCG", "UAGC") # Translation table
                # Create each of the 6 reading frames
                self.reading_frame_1 = "".join([self.amino_acid_code[self.standard_genetic_code[codon]] 
                                                for codon in re.findall(".{3}", self.RNA)])
                self.reading_frame_2 = "".join([self.amino_acid_code[self.standard_genetic_code[codon]] 
                                                for codon in re.findall(".{3}", self.RNA[1:-2])])
                self.reading_frame_3 = "".join([self.amino_acid_code[self.standard_genetic_code[codon]] 
                                                for codon in re.findall(".{3}", self.RNA[2:-1])])
                self.reading_frame_4 = "".join([self.amino_acid_code[self.standard_genetic_code[codon]] 
                                                for codon in re.findall(".{3}", 
                                                                        self.RNA.translate(self.translation_table))[::-1]])
                self.reading_frame_5 = "".join([self.amino_acid_code[self.standard_genetic_code[codon]] 
                                                for codon in re.findall(".{3}", 
                                                                        self.RNA[2:-1].translate(self.translation_table))[::-1]])
                self.reading_frame_6 = "".join([self.amino_acid_code[self.standard_genetic_code[codon]] 
                                                for codon in re.findall(".{3}", 
                                                                        self.RNA[1:-2].translate(self.translation_table))[::-1]])
                # Put all the frames into a list
                self.all_reading_frames = [self.reading_frame_1, 
                                           self.reading_frame_2, 
                                           self.reading_frame_3, 
                                           self.reading_frame_4, 
                                           self.reading_frame_5, 
                                           self.reading_frame_6]
                # Create a tuple for each reading frame of tuples for each start and stop codons match
                self.start_and_stop_codons = tuple(tuple("M" + match + "_" for match in re.findall("M(.*?)_", reading_frame)) for reading_frame in self.all_reading_frames)
            # If there are non-nucleotides, halt
            else:
                raise Exception("The sequence is not valid because in {0}, the number of non-nucleotides is {1}.".format(filename, self.nucleotide_totals["N"]))
    
    # Method for translation into selected reading frames
    def translateDNA(self):
        keep_going = True
        while keep_going: # Exit the loop based on user input
            print("You will be prompted to enter a specific reading frame.")
            print("Valid frames are 1, 2, 3, 4, 5, and 6.")
            print("If you would like to see all reading frames, you may enter 0.")
            print("="*20 + " WARNING: ALL OTHER INPUTS WILL CAUSE THE PROGRAM TO HALT " + "="*20)
            user_input = input("Enter a specific reading frame: ")
            if user_input == "0":
                print("\n\n".join(self.all_reading_frames))
                return self.all_reading_frames
            elif user_input == "1":
                print(self.reading_frame_1)
                return self.reading_frame_1
            elif user_input == "2":
                print(self.reading_frame_2)
                return self.reading_frame_2
            elif user_input == "3":
                print(self.reading_frame_3)
                return self.reading_frame_3
            elif user_input == "4":
                print(self.reading_frame_4)
                return self.reading_frame_4
            elif user_input == "5":
                print(self.reading_frame_5)
                return self.reading_frame_5
            elif user_input == "6":
                print(self.reading_frame_6)
                return self.reading_frame_6
            else:
                # None of the valid options were used
                print("INVALID INPUT; THE PROGRAM WILL NOW HALT.")
                keep_going = False
                
    def predictRF(self):
        # The tuple of tuples was already created in the constructor so just return it
        return self.start_and_stop_codons
        
    def longestRF(self):
        # Use a generator comprehension to find the indices of the longest start and stop codon matches
        max_value, (max_outer, max_inner) = max(((start_and_stop, (i, j))
                                                 for i, reading_frame in enumerate(self.start_and_stop_codons)
                                                 for j, start_and_stop in enumerate(reading_frame)), 
                                                key=lambda x: max(len(e) for e in x))
        return max_outer + 1, max_value # Return a tuple with two elements: the reading frame and the longest match

# Set up main function for tests
def main():
    # Load the example files
    sequence_1 = Sequence("data/RFdata.txt")
    sequence_2 = Sequence("data/RFdata2.txt")
    sequence_3 = Sequence("data/RFdata3.txt")

    # Show the statistics for the first file
    print("Statistics for " + sequence_1.filename + ":")
    print(sequence_1.nucleotide_totals)

    # Show the longest longest start and stop codons match 
    # and its associated reading frame for the first file
    print("Example use of longestRF method on " + sequence_1.filename + ":")
    print(sequence_1.longestRF())

    # Show the statistics for the second file
    print("Statistics for " + sequence_2.filename + ":")
    print(sequence_2.nucleotide_totals)

    # Show the longest longest start and stop codons match 
    # and its associated reading frame for the second file
    print("Example use of longestRF method on " + sequence_2.filename + ":")
    print(sequence_2.longestRF())

    # Show the statistics for the third file
    print("Statistics for " + sequence_3.filename + ":")
    print(sequence_3.nucleotide_totals)

    # Show the longest longest start and stop codons match 
    # and its associated reading frame for the third file
    print("Example use of longestRF method on " + sequence_3.filename + ":")
    print(sequence_3.longestRF())

    # Test translation based on user input
    print("Example use of translateDNA method on " + sequence_3.filename + ":")
    result = sequence_3.translateDNA()
    
    # Show all potential start and stop codon matches for all reading frames
    print("Example use of predictRF method on " + sequence_3.filename + ":")
    print(sequence_3.predictRF())

# Run the main function along with all the tests
main()
