### CMSC220 Assignment 4
The fourth assignment for CMSC 220: Bioinformatics, which involves using Python to identify the origin of replication in DNA.

#### Program Output:
```
K-Mers for Vibrio Cholerae, with reading frames:
------------------------------------------------
1|3|ATC         GAT         TCT         TTG
 |6|GATCAA      ATCAAT      CGTAAG      CTTCTA
 |9|ATCAATGAT   CAACGTAAG   CTTCTAAGC   ATGATCAAG
 ================================================
2|3|TGA         CTC         CTT         TCA
 |6|TGACAT      CTCTTG      TCAATG      ATCAAC
 |9|TCAATGATC   AACGTAAGC   TTCTAAGCA   TGATCAAGG
 ================================================
3|3|TGA         GAT         ATC         TTT
 |6|GACATC      TACTCT      TGGCCA      AATGAA
 |9|CAATGATCA   ACGTAAGCT   TCTAAGCAT   GATCAAGGT
 ================================================
4|3|TAG         CTA         AGA         AAC
 |6|CTAGTT      GCAAAG      CTAGTA      CGAGAA
 |9|GAACTAGTA   CGACGACGA      CTAGTT   CCTTCTTA
 ================================================
5|3|ACT         AAA         CTA         TAG
 |6|TTACTT      ATGAGA      ACCGGT      CTGTAG
 |9|ACTAGTAGC   ACGACGAGA    CTAGTTCG   TTCTTA
 ================================================
6|3|ACT         GAG         GAA         TAG
 |6|GAGAAC      ACTGTA      TAGTAG      GACGAC
 |9|AACTAGTAG   GACGACGAG     CTAGTTC   CTTCTTA
 ================================================

K-Mers for Vibrio Cholerae, without reading frames:
---------------------------------------------------
1|3|TGA         ATC         GAT         TCA
 |6|TGATCA      ATGATC      GATCAA      ATCAAG
 |9|CTCTTGATC   TCTTGATCA   CTTGATCAT   AAGCATGAT
 ================================================

K-Mers for Thermotoga Petrophila, with reading frames:
------------------------------------------------------
1|3|AAA         ATT         TAC         ATA
 |6|TGGTAG      AAACAA      CTTACC      AACTCT
 |9|AACTCTATA   CCTCCTTTT   TGTCGAATT   TGTGTGATT
 ================================================
2|3|TTT         ACC         TTA         GAA
 |6|GGTAGG      ACTCTA      TACCTC      CTTTTT
 |9|ACTCTATAC   CTCCTTTTT   GTCGAATTT   GTGTGATTT
 ================================================
3|3|AAA         TTG         AAT         TAT
 |6|GTAGGT      TTCAGA      TACCAC      CTCTAT
 |9|CTCTATACC   TCCTTTTTG   TCGAATTTG   TGTGATTTA
 ================================================
4|3|TTT         TAA         ATG         TAT
 |6|GAATGG      TTTGTT      ACCATC      CACCAT
 |9|CTCCACCAT   TAGACTTTT   CATATTAAC   TATTATCGT
 ================================================
5|3|TTT         TTA         AAC         TGA
 |6|ATGGTG      AAGTCT      CATCCA      CCATTT
 |9|CCACCATTT   GACTTTTCT   TATTAACTA   TTATCGTCA
 ================================================
6|3|AAA         TGG         AAT         CTT
 |6|CCATCC      ACCATT      TTCTCC      AGACTT
 |9|TCCACCATT   AGACTTTTC   ATATTAACT   ATTATCGTC
 ================================================

K-Mers for Thermotoga Petrophila, without reading frames:
---------------------------------------------------------
1|3|AAA         TTT         ATT         ACC
 |6|TACCAC      ACCTAC      CCTACC      CTACCA
 |9|ACCTACCAC   AAACCTACC   AACCTACCA   CCTACCACC
 ================================================
```