# Noah Segal-Gould <ns2349@bard.edu>
# 2018-02-13
# CMSC 220
# Lab Assignment 2
# My code is based in part on work done collectively in class as well as code provided on the course Moodle 2.

import re

def main(protein_filename="protein.txt", 
         nucleotide_filename="nucleotide.txt", 
         pdb_filename="protein.pdb"):
    
    protein_text = open(protein_filename).read()
    nucleotide_text = open(nucleotide_filename).read()
    pdb_text = open(pdb_filename).read()
        
    print("Protein Name: {}\n".format(re.search(r"\[protein=(.*?)\]", protein_text).group(1)))
    
    keep_going = True
    
    print("Please enter one of the following options:\n")
    print("1: Nucleotide data:\n\tThis option will output the number of each nucleotide in the sequence and the frequency of each nucleotide.")
    print("2: Protein data:\n\tThis option will allow you to input a specific motif to find, and display the location of each motif in the data.")
    print("3: Centroid data:\n\tThis option will output the centroid location of the protein.\n")
    print("="*20 + " WARNING " + "="*20)
    print("If you fail to input one of these options, the program will halt.")
        
    while(keep_going):
        
        user_input = input("\nPlease enter either 1, 2, or 3: ")
        
        if user_input == "1":
            
            a_regex = re.compile(r"[Aa]")
            t_regex = re.compile(r"[Tt]")
            c_regex = re.compile(r"[Cc]")
            g_regex = re.compile(r"[Gg]")
            non_nucleotide_regex = re.compile(r"[^ATCGatcg]")
            
            nucleotides = "".join([line for line in nucleotide_text.split("\n") if ">" not in line])
            nucleotides_length = len(nucleotides)
            
            number_of_adenine = len(a_regex.findall(nucleotides))
            frequency_of_adenine = "{0:.0f}%".format(number_of_adenine/nucleotides_length * 100)
            
            number_of_thymine = len(t_regex.findall(nucleotides))
            frequency_of_thymine = "{0:.0f}%".format(number_of_thymine/nucleotides_length * 100)
            
            number_of_cytosine = len(c_regex.findall(nucleotides))
            frequency_of_cytosine = "{0:.0f}%".format(number_of_cytosine/nucleotides_length * 100)
            
            number_of_guanine = len(g_regex.findall(nucleotides))
            frequency_of_guanine = "{0:.0f}%".format(number_of_guanine/nucleotides_length * 100)
            
            number_of_non_nucleotides = len(non_nucleotide_regex.findall(nucleotides))
            
            output_list = [number_of_adenine, frequency_of_adenine, 
                           number_of_thymine, frequency_of_thymine, 
                           number_of_cytosine, frequency_of_cytosine, 
                           number_of_guanine, frequency_of_guanine]
            
            print("A: {0} ({1}), T: {2} ({3}), C: {4} ({5}), G: {6} ({7})".format(*output_list))
            print("Total number of unidentified nucleotides: {}".format(number_of_non_nucleotides))
            
        elif user_input == "2":
            
            protein = "".join([line for line in protein_text.split("\n") if ">" not in line])

            motif_user_input = input("Please enter a motif: ").upper()
            motif_locations_list = [match.start() for match in re.finditer(motif_user_input, protein)]
            
            if motif_locations_list:
                print("The following locations were found for that motif (\"{0}\"):\n\t{1}".format(motif_user_input, 
                                                                                                   str(motif_locations_list)[1:-1]))
            else:
                print("\"{}\" not found; returning to initial prompt.".format(motif_user_input))
        
        elif user_input == "3":
            
            pdb = pdb_text.split("\n")            
            natoms = 0
            
            for line in pdb:
                if line[:6] == "ATOM  ":
                    natoms += 1
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
            if natoms != 0:
                centroid_position = str([round(x/natoms, 6), 
                                         round(y/natoms, 6), 
                                         round(z/natoms, 6)])[1:-1]
                print("Centroid size: {0}\nCentroid position: {1}".format(natoms, centroid_position))                
            else:
                print("No centroid found; returning to initial prompt.")
        
        else:
            print("Invalid input; this program will now halt.")
            keep_going = False

main()