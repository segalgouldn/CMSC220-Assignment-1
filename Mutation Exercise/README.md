### CMSC220 Mutation Exercise
An in-class exercise for CMSC 220: Bioinformatics, which involves using Python to create and identify nucleotide sequence mutations.

#### Program Output:
```
==================== Testing Mutations ====================
The sequences are ATCGATCGAT & TACGGACTCA
The sequence "ATCGATCGAT" undergoes a point mutation and becomes "ATCGATTGAT"
The sequence "ATCGATCGAT" undergoes a deletion mutation and becomes "ATCGACGAT"
The sequence "ATCGATCGAT" undergoes an insertion mutation and becomes "ATCGATCCGAT"
The sequence "TACGGACTCA" undergoes a point mutation and becomes "TATGGACTCA"
The sequence "TACGGACTCA" undergoes a deletion mutation and becomes "TAGGACTCA"
The sequence "TACGGACTCA" undergoes an insertion mutation and becomes "TACGGGACTCA"
==================== Testing Mutation Detection ====================
Going from ATCGATCGAT to ATCGATTGAT there was a missense mutation
Going from ATCGATCGAT to ATCGACGAT there was a frameshift deletion mutation
Going from ATCGATCGAT to ATCGATCCGAT there was a frameshift insertion mutation
Going from TACGGACTCA to TATGGACTCA there was a silent mutation
Going from TACGGACTCA to TAGGACTCA there was a frameshift deletion mutation
Going from TACGGACTCA to TACGGGACTCA there was a frameshift insertion mutation
```