# Noah Segal-Gould <ns2349@bard.edu>
# 2018-02-13
# CMSC 220
# Lab Assignment 2
# My code is based in part on work done collectively in class as well as code provided on the course Moodle 2.

# Import Python's regular expressions library
import re

# Define a main function which will be run at the end of the file
def main(protein_filename="protein.txt", # Make a default protein filename 
         nucleotide_filename="nucleotide.txt", # Make a default nucleotide filename
         pdb_filename="protein.pdb"): # Make a default PDB filename
    
    # Open the three files and get the text inside as one large string for each
    protein_text = open(protein_filename).read()
    nucleotide_text = open(nucleotide_filename).read()
    pdb_text = open(pdb_filename).read()
    
    # Use a regular expression to find the name of the protein
    print("Protein Name: {}\n".format(re.search(r"\[protein=(.*?)\]", protein_text).group(1)))
    
    # Prepare the program's main while loop to run until a variable changes
    keep_going = True
    
    # Print instructions
    print("Please enter one of the following options:\n")
    print("1: Nucleotide data:\n\tThis option will output the number of each nucleotide in the sequence and the frequency of each nucleotide.")
    print("2: Protein data:\n\tThis option will allow you to input a specific motif to find, and display the location of each motif in the data.")
    print("3: Centroid data:\n\tThis option will output the centroid location of the protein.\n")
    print("="*20 + " WARNING " + "="*20)
    print("If you fail to input one of these options, the program will halt.")
    
    # The main while loop; it executes until the user inputs something other than the three main options
    while(keep_going):
        
        # Take user input as a string
        user_input = input("\nPlease enter either 1, 2, or 3: ")
        
        # Check if the first option was used
        if user_input == "1":
            
            # Prepare regular expression patterns upon which to use findall, for use in finding some statistics
            a_regex = re.compile(r"[Aa]")
            t_regex = re.compile(r"[Tt]")
            c_regex = re.compile(r"[Cc]")
            g_regex = re.compile(r"[Gg]")
            
            # Prepare regular expression pattern for identification of errors
            non_nucleotide_regex = re.compile(r"[^ATCGatcg]")
            
            # Combine together all nucleotide lines from the file
            nucleotides = "".join([line for line in nucleotide_text.split("\n") if ">" not in line])
            nucleotides_length = len(nucleotides)
            
            # Calculate nucleotide counts and frequencies for A, T, C, and G
            number_of_adenine = len(a_regex.findall(nucleotides))
            frequency_of_adenine = "{0:.0f}%".format(number_of_adenine/nucleotides_length * 100)
            
            number_of_thymine = len(t_regex.findall(nucleotides))
            frequency_of_thymine = "{0:.0f}%".format(number_of_thymine/nucleotides_length * 100)
            
            number_of_cytosine = len(c_regex.findall(nucleotides))
            frequency_of_cytosine = "{0:.0f}%".format(number_of_cytosine/nucleotides_length * 100)
            
            number_of_guanine = len(g_regex.findall(nucleotides))
            frequency_of_guanine = "{0:.0f}%".format(number_of_guanine/nucleotides_length * 100)
            
            # Calculate number of errors
            number_of_non_nucleotides = len(non_nucleotide_regex.findall(nucleotides))
            
            output_list = [number_of_adenine, frequency_of_adenine, 
                           number_of_thymine, frequency_of_thymine, 
                           number_of_cytosine, frequency_of_cytosine, 
                           number_of_guanine, frequency_of_guanine]
            
            # Print the statistics
            print("A: {0} ({1}), T: {2} ({3}), C: {4} ({5}), G: {6} ({7})".format(*output_list))
            print("Total number of unidentified nucleotides: {}".format(number_of_non_nucleotides))
        
        # Check if the second option was used
        elif user_input == "2":
            
            # Same as the one for nucleotides; joins together all lines that do not contain ">"
            protein = "".join([line for line in protein_text.split("\n") if ">" not in line])
                  
            # Take user input and capitalize it
            motif_user_input = input("Please enter a motif: ").upper()
            # Use regular expressions to find indices of substrings
            motif_locations_list = [match.start() for match in re.finditer(motif_user_input, protein)]
            
            # Print the results
            if motif_locations_list:
                print("The following locations were found for that motif (\"{0}\"):\n\t{1}".format(motif_user_input, 
                                                                                            str(motif_locations_list)[1:-1])) # Convert the list to a comma-separated string
            # Print the error
            else:
                print("\"{}\" not found; returning to initial prompt.".format(motif_user_input))
        
        # Check if the third and final option was used
        elif user_input == "3":
            
            # What follows in this case is mostly adapted from the centroid code shown in class and posted on Moodle 2
            pdb = pdb_text.split("\n")            
            natoms = 0
            
            for line in pdb:
                if line[:6] == "ATOM  ":
                    natoms += 1
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    xsum += x
                    ysum += y
                    zsum += z
                    
            if natoms != 0:
                centroid_position = str([round(xsum/natoms, 6), 
                                         round(ysum/natoms, 6), 
                                         round(zsum/natoms, 6)])[1:-1]
                print("Centroid size: {0}\nCentroid position: {1}".format(natoms, centroid_position))                
            else:
                print("No centroid found; returning to initial prompt.")
        
        # The alternative to the last option is that the program gives up
        else:
            print("Invalid input; this program will now halt.")
            # The while loop should exit since the user didn't input anything meaningful
            keep_going = False

# Call the main function and run the program
main()
