### CMSC220 Assignment 4
The fourth assignment for CMSC 220: Bioinformatics, which involves using Python to identify the origin of replication in DNA.

#### Program Output:
```
K-Mers for Vibrio Cholerae, with reading frames:
------------------------------------------------
1|3|TGA         GAT         CTT         CTC
 |6|GATCAA      TGGCCA      ATCAAT      CGTAAG
 |9|ATCAATGAT   CAACGTAAG   CTTCTAAGC   ATGATCAAG
 ================================================
2|3|ATC         GAT         CAT         TCT
 |6|ATCAAG      ACTTGT      CATGAT      TCAATG
 |9|TCAATGATC   AACGTAAGC   TTCTAAGCA   TGATCAAGG
 ================================================
3|3|TGA         TCA         ATG         CAA
 |6|TGATCA      CAATGA      TCAACG      TAAGCT
 |9|CAATGATCA   ACGTAAGCT   TCTAAGCAT   GATCAAGGT
 ================================================
4|3|ACT         CTA         GAA         GAG
 |6|CTAGTT      ACCGGT      GCAAAG      CTAGTA
 |9|GTAGCAAAG   CGAGAACTA   GTTCGACGA   TTCTTACTA
 ================================================
5|3|AGT         ACT         TAC         AAA
 |6|ACTAGT      AGTAGC      AGAACT      ACGACG
 |9|AGAACTAGT   TCGACGACG   CTTACTAGT   AAATGCCTT
 ================================================
6|3|TAG         CTA         GTA         AGA
 |6|TAGTTC      GTACTA      TGAACA      TAGTAG
 |9|GAGAACTAG   TTCGACGAC   TCTTACTAG   AAAATGCCT
 ================================================

K-Mers for Vibrio Cholerae, without reading frames:
---------------------------------------------------
1|3|TGA         ATC         GAT         TCA
 |6|TGATCA      ATGATC      GATCAA      ATCAAG
 |9|ATGATCAAG   CTCTTGATC   TCTTGATCA   CTTGATCAT
 ================================================

K-Mers for Thermotoga Petrophila, with reading frames:
------------------------------------------------------
1|3|AAC         TTT         AAA         TAT
 |6|AGGTTT      GGTGGT      CCTACC      GTGGTA
 |9|AACTCTATA   CCTCCTTTT   TGTCGAATT   TGTGTGATT
 ================================================
2|3|ATT         TAC         ACT         TTT
 |6|GTGGTA      AATTGA      ACTTAC      CTACCA
 |9|ACTCTATAC   CTCCTTTTT   GTCGAATTT   GTGTGATTT
 ================================================
3|3|AAA         ACC         ATT         TTT
 |6|TGGTAG      CTGAAA      CTTACC      ATTTCA
 |9|CTCTATACC   TCCTTTTTG   TCGAATTTG   TGTGATTTA
 ================================================
4|3|TTT         AAA         TTG         CAT
 |6|TCCAAA      CACCAT      GGATGG      CCACCA
 |9|TTTCTCCAC   AACTAGACT   CGTCATATT   GATTATTAT
 ================================================
5|3|TTT         TGG         TAA         AAA
 |6|ACCATC      GACTTT      ATGGTG      TAAAGT
 |9|TCTCCACCA   CTAGACTTT   TCATATTAA   TTATTATCG
 ================================================
6|3|TAA         ATG         TGA         AAA
 |6|ACCATT      TTAACT      GATGGT      TGAATG
 |9|TTCTCCACC   ACTAGACTT   GTCATATTA   ATTATTATC
 ================================================

K-Mers for Thermotoga Petrophila, without reading frames:
---------------------------------------------------------
1|3|AAA         TTT         ATT         ACC
 |6|TACCAC      ACCTAC      CCTACC      CTACCA
 |9|ACCTACCAC   TGGTAGGTT   GGTAGGTTT   AAACCTACC
 ================================================
```