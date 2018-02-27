### CMSC220 Assignment 2
The second assignment for CMSC 220: Bioinformatics, which involves using Python to find statistics about nucleotides and proteins using input files.

#### Program Output:
```
Protein Name: cytochrome b

Please enter one of the following options:

1: Nucleotide data:
	This option will output the number of each nucleotide in the sequence and the frequency of each nucleotide.
2: Protein data:
	This option will allow you to input a specific motif to find, and display the location of each motif in the data.
3: Centroid data:
	This option will output the centroid location of the protein.

==================== WARNING ====================
If you fail to input one of these options, the program will halt.

Please enter either 1, 2, or 3: 1
A: 1839 (20%), T: 3250 (35%), C: 2050 (22%), G: 2122 (23%)
Total number of unidentified nucleotides: 0

Please enter either 1, 2, or 3: 2
Please enter a motif: LAV
The following locations were found for that motif ("LAV"):
	573, 576, 1563, 2126, 2318, 2525

Please enter either 1, 2, or 3: 3
Centroid size: 37804
Centroid position: 23.065753, -0.07208, 22.827757

Please enter either 1, 2, or 3: exit
Invalid input; this program will now halt.
```
