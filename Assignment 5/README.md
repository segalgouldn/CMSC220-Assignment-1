### CMSC220 Assignment 5
The fifth assignment for CMSC 220: Bioinformatics, which involves using Python to cause and detect genomic variations.

#### Program Output:
```
* Testing SequenceVariation.multiplePointMutations

seq1: GCACGTATTGATTGGCCTGTACCTA
seq2: TCACGAATTCAGAGTCCAGTACCCG

seq1 and seq2 have 9 differences

* Testing SequenceVariation.doCrossover

seq1: GCACGTATTGATTGGCCTGTACCTA and seq2: TCACGAATTCAGAGTCCAGTACCCG
become
seq3: GCACGTATTGATTGGCCTGTACCCG and seq4: TCACGAATTCAGAGTCCAGTACCTA

* Testing SequenceVariation.doCopyNumberVariationAddition

Copy-number variation addition on
GCACGTATTGATTGGCCTGTACCTA
produces
GCACGTATTGATTGATTGGCCTGTACCTA

* Testing SequenceVariation.doCopyNumberVariationSubtraction

Copy-number variation subtraction on
GCACGTATTGATTGGCCTGTACCTA
produces
GCACGTATTGGCCTGTACCTA

* Testing SequenceVariation.doTransposon

Transposon variation on
GCACGTATTGATTGGCCTGTACCTA
produces
GCACGTGGTTGCACGTATTGATTGGCCTGTACCTA

* Testing detection of mutations:

Going from GCACGTATTGATTGGCCTGTACCTA to GCACGTATTGATTGGCCTGTACCCG there was a crossover mutation
Going from GCACGTATTGATTGGCCTGTACCTA to TCACGAATTCAGAGTCCAGTACCTA there was a crossover mutation
Going from GCACGTATTGATTGGCCTGTACCTA to GCACGTATTGATTGATTGGCCTGTACCTA there was a copy variation mutation
Going from GCACGTATTGATTGGCCTGTACCTA to GCACGTATTGGCCTGTACCTA there was a copy variation mutation
Going from GCACGTATTGATTGGCCTGTACCTA to GCACGTGGTTGCACGTATTGATTGGCCTGTACCTA there was a transposon insertion mutation

* Testing similarity scores:

Similarity between GCACGTATTGATTGGCCTGTACCTA and TCACGAATTCAGAGTCCAGTACCCG: -11.0
Similarity between GCACGTATTGATTGGCCTGTACCTA and GCACGTATTGATTGGCCTGTACCCG: 17.0
Similarity between GCACGTATTGATTGGCCTGTACCTA and TCACGAATTCAGAGTCCAGTACCTA: -3.0
Similarity between GCACGTATTGATTGGCCTGTACCTA and GCACGTATTGATTGATTGGCCTGTACCTA: -7.0
Similarity between GCACGTATTGATTGGCCTGTACCTA and GCACGTATTGGCCTGTACCTA: -11.0
Similarity between GCACGTATTGATTGGCCTGTACCTA and GCACGTGGTTGCACGTATTGATTGGCCTGTACCTA: -35.0

* Testing multiple sequence alignment:

0 -GCACGT-ATTG-A--T-T-G-G---CCTGTACCTA
1 T-CACGA-ATTC-A--G-A-G-T---CCAGTACCCG
2 -GCACGT-ATTG-A--T-T-G-G---CCTGTACCCG
3 T-CACGA-ATTC-A--G-A-G-T---CCAGTACCTA
4 -GCACGT-ATTG-A--T-T-GATTGGCCTGTACCTA
5 -----------GCA-CG-T--ATTGGCCTGTACCTA
6 -GCACGTGGTTGCACGTATTGATTGGCCTGTACCTA
```