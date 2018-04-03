# Noah Segal-Gould, Aaron Krapf, and Zoe Terhune
# 2018-04-03
# The three of us collaborated on this assignment.
# We made use of the textbook as well as code on Moodle 2.

from random import randint, choice

# From Alignments.py on Moodle 2
DNA_2 = {'G': { 'G': 1, 'C':-3, 'A':-3, 'T':-3, 'N':0 },
         'C': { 'G':-3, 'C': 1, 'A':-3, 'T':-3, 'N':0 },
         'A': { 'G':-3, 'C':-3, 'A': 1, 'T':-3, 'N':0 },
         'T': { 'G':-3, 'C':-3, 'A':-3, 'T': 1, 'N':0 },
         'N': { 'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N':0 }}

# From consensusMultAlignClass.py on Moodle 2
def sequenceAlign(seqA, seqB, simMatrix=DNA_2, insert=8, extend=4):
  
  numI = len(seqA) + 1
  numJ = len(seqB) + 1
  
  scoreMatrix = [[0] * numJ for x in range(numI)]
  routeMatrix = [[0] * numJ for x in range(numI)]
  
  for i in range(1, numI):
    routeMatrix[i][0] = 1
  
  for j in range(1, numJ):
    routeMatrix[0][j] = 2
  
  for i in range(1, numI):
    for j in range(1, numJ):
      
      penalty1 = insert
      penalty2 = insert
      
      if routeMatrix[i-1][j] == 1:
        penalty1 = extend
      
      elif routeMatrix[i][j-1] == 2:
        penalty2 = extend
      
      similarity = simMatrix[ seqA[i-1] ][ seqB[j-1] ]
      
      paths = [scoreMatrix[i-1][j-1] + similarity, # Route 0
               scoreMatrix[i-1][j] - penalty1, # Route 1
               scoreMatrix[i][j-1] - penalty2] # Route 2
      
      best = max(paths)
      route = paths.index(best)
      
      scoreMatrix[i][j] = best
      routeMatrix[i][j] = route
  
  alignA = []
  alignB = []
  
  i = numI-1
  j = numJ-1
  score = scoreMatrix[i][j]
    
  
  while i > 0 or j > 0:
    route = routeMatrix[i][j]
    
    if route == 0: # Diagonal
      alignA.append( seqA[i-1] )
      alignB.append( seqB[j-1] )
      i -= 1
      j -= 1
    
    elif route == 1: # Gap in seqB
      alignA.append( seqA[i-1] )
      alignB.append( '-' )
      i -= 1
    
    elif route == 2: # Gap in seqA
      alignA.append( '-' )
      alignB.append( seqB[j-1] )
      j -= 1
  
  alignA.reverse()
  alignB.reverse()
  alignA = ''.join(alignA)
  alignB = ''.join(alignB)
  
  return score, alignA, alignB

# From consensusMultAlignClass.py on Moodle 2
def consensus(alignment, threshold=0.25):
  
  #finding length of first sequence
  n = len(alignment[0])
  #number of sequences
  nSeq = float(len(alignment))
  consensus = ''
  
  #loop through all of the bases in the sequence
  for i in range(n):
    #creates new dictionary
    counts = {}
    
    #for loop for every sequence in alignment
    for seq in alignment:
      #sets letter to ith nucleotide in sequence
      letter = seq[i]
      #if nucleotide is dash (gap), skip it
      if letter == '-':
        continue
    #sets value of letter to one plus previous
    #if previous doesn't exist, makes new key and sets its value to 1
      counts[letter] = counts.get(letter, 0) + 1
    
    #new list
    fractions = []
    #for every key in dictionary, creates a fraction
    for letter in counts:
        #take value/number of sequences
      frac = counts[letter]/nSeq
      #appends a list of lists with the fraction and letter
      #frac is the average/percetage of each letter in the sequences
      fractions.append([frac, letter])
    
    #sorts the list in order from leat frequent to most
    fractions.sort()
    #pull out last entry (most frequent) and set for bestFrac, bestLet
    bestFraction, bestLetter = fractions[-1]
    
    #testing if frantion is less than or equal to threshold
    #if so add 'N' - no consensus
    if bestFraction <=  threshold:
      consensus += 'N'
    #otherwise adds bestletter to consensus
    else:
      consensus += bestLetter
  
  return consensus

# From consensusMultAlignClass.py on Moodle 2
#inputs are the seqs, the threshold, and a similarity Matrix
def consensusMultipleAlign(seqs, threshold, simMatrix):
  #counts number of sequences
  n = len(seqs)
  #setting default penalty and new list
  insert = 2;
  extend = 4;
  multipleAlign = []
  
  i = 0
  #looping through sequences
  #set j as i+1 to compare first to second seq
  for j in range(i+1,n):
    #seqB starts as first sequence
    seqB = seqs[j]
    
    #at the beginning, seqA is our first sequence
    if not multipleAlign:
      seqA = seqs[i]
      #Align our two sequences (with gaps!)
      #outputs: score and two newly aligned sequences
      score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix,insert, extend)
      #adds newly aligned sequences to our list
      multipleAlign.append(alignA)
      multipleAlign.append(alignB)
    
    
    else:
      # seqA is our consensus sequence and align seqB to our consensus seq
      seqA = consensus(multipleAlign, threshold)
      score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix,insert, extend)
      #creates new list
      gaps = []
      #for everythin in alignA list, k is index, and letter is value
      for k, letter in enumerate(alignA):
        # if letter is a dash (there's a gap)
        if letter == '-':
          #appends index to gap
          #keeps track of locations of gaps
          gaps.append(k)
      
      #enumerate through each sequence in multipleAlign
      for k, seq in enumerate(multipleAlign):
        #go through list of gaps
        for gap in gaps:
          # add a dash at each index for each gap
          seq = seq[:gap] + '-' + seq[gap:]
        
        #
        multipleAlign[k] = seq
      
      #
      multipleAlign.append(alignB)
  
  
  #
  for k, seq in enumerate(multipleAlign):
    print(k, seq)

# Adapted from the textbook
def similarity(sequence_1, sequence_2):
    numPlaces = min(len(sequence_1), len(sequence_2))
    totalScore = 0.0
    for i in range(numPlaces):
        residueA = sequence_1[i]
        residueB = sequence_2[i]
        totalScore += DNA_2[residueA][residueB]
    return totalScore

# To be used for identifying mutations
def sequence_difference(sequence_1, sequence_2):
    if len(sequence_1) == len(sequence_2):
        differences = 0
        for i in range(len(sequence_1)):
            if sequence_1[i] != sequence_2[i]:
                differences += 1
        return differences
    else:
        raise Exception("The lengths of {0} and {1} are not equal.".format(sequence_1, sequence_2))

class SequenceVariation():
    def __init__(self, sequence="GCACGTATTGATTGGCCTGTACCTA"):
        self.dna = "".join([n if n in "ATCG" else "N" for n in sequence.upper()])
    
    def multiplePointMutations(self, n=10):
        def doPointMutation(sequence):
            default_nucleotides = "ATCG"
            mutation_index = randint(0, len(sequence) - 1)
            mutation_index_nucleotide = sequence[mutation_index]
            
            mutation_beginning = sequence[:mutation_index]
            mutation_middle = choice(default_nucleotides.replace(mutation_index_nucleotide, ""))
            mutation_end = sequence[mutation_index + 1:]
            
            mutation_final = mutation_beginning + mutation_middle + mutation_end
            
            if mutation_final == sequence:
                mutation_final = doPointMutation(sequence)
            
            return mutation_final
        
        if n > len(self.dna) or not isinstance(n, int):
            raise Exception("Invalid Mutation Number: {}".format(n))
        
        newSeq = self.dna
        for i in range(n):
            newSeq = doPointMutation(newSeq)
        return newSeq
    
    def doCrossover(self, sequence):
        sequence_length = len(sequence)
        dna_length = len(self.dna)
        if sequence_length != dna_length:
            raise Exception("Can not perform crossover on strings of different lengths.")
        
        crossover_index = randint(1, sequence_length - 2)
        
        sequence_first, sequence_second = sequence[:crossover_index], sequence[crossover_index:]
        dna_first, dna_second = self.dna[:crossover_index], self.dna[crossover_index:]
        
        return dna_first + sequence_second, sequence_first + dna_second
    
    def doCopyNumberVariationAddition(self):
        for i in range(3, (len(self.dna)//2) + 1):
            kmers = [self.dna[j:j + i] for j in range(len(self.dna) - i + 1)]
            for k, el in enumerate(kmers):
                if k + i < len(kmers) and el == kmers[k + i]:
                    match_index = k + (i * 2)
                    return self.dna[:match_index] + el + self.dna[match_index:]
    
    def doCopyNumberVariationSubtraction(self):
        for i in range(3, (len(self.dna)//2) + 1):
            kmers = [self.dna[j:j + i] for j in range(len(self.dna) - i + 1)]
            for k, el in enumerate(kmers):
                if k + i < len(kmers) and el == kmers[k + i]:
                    match_index = k + (i * 2)
                    return self.dna[:match_index - i] + self.dna[match_index:]
    
    def doTransposon(self, transposonSequence="ACGTGGTTGCACGT"):
        return self.dna.replace(transposonSequence[:4], transposonSequence, 2)

class SequenceVariationDetection():
    def __init__(self, sequence_1, sequence_2):
        self.sequence_1 = sequence_1
        self.sequence_2 = sequence_2
    
    def detect_variation(self):
        sequence_1 = self.sequence_1
        sequence_2 = self.sequence_2
        if sequence_1 != sequence_2:
            if len(sequence_1) == len(sequence_2) and sequence_difference(sequence_1, sequence_2) > 1:
                return "Going from {0} to {1} there was a crossover mutation".format(sequence_1, sequence_2)
            elif len(sequence_1) == len(sequence_2):
                return "Going from {0} to {1} there was a point mutation".format(sequence_1, sequence_2)
            elif abs(len(sequence_1) - len(sequence_2)) == len("ACGTGGTTGCACGT") - 4:
                return "Going from {0} to {1} there was a transposon insertion mutation".format(sequence_1, sequence_2)
            elif abs(len(sequence_1) - len(sequence_2)) == 4:
                return "Going from {0} to {1} there was a copy variation mutation".format(sequence_1, sequence_2)
            elif len(sequence_2) > len(sequence_1):
                return "Going from {0} to {1} there was a frameshift insertion mutation".format(sequence_1, sequence_2)
            elif len(sequence_2) < len(sequence_1):
                return "Going from {0} to {1} there was a frameshift deletion mutation".format(sequence_1, sequence_2)
            else:
                raise Exception("No mutation occurred going from {0} to {1}".format(sequence_1, sequence_2))
        else:
            raise Exception("No mutation occurred going from {0} to {1}".format(sequence_1, sequence_2))

#Tests:

sequenceVariationObject = SequenceVariation()
print("* Testing SequenceVariation.multiplePointMutations\n")

seq1 = "GCACGTATTGATTGGCCTGTACCTA"
print("seq1: {}".format(seq1))

seq2 = SequenceVariation().multiplePointMutations()
print("seq2: {}".format(seq2))

print("\nseq1 and seq2 have {} differences".format(sequence_difference(seq1, seq2)))

print("\n* Testing SequenceVariation.doCrossover\n")

seq3, seq4 = sequenceVariationObject.doCrossover(seq2)

print("seq1: {0} and seq2: {1}\nbecome\nseq3: {2} and seq4: {3}\n".format(seq1, seq2, seq3, seq4))

print("* Testing SequenceVariation.doCopyNumberVariationAddition\n")

seq5 = sequenceVariationObject.doCopyNumberVariationAddition()

print("Copy-number variation addition on \n{0}\nproduces\n{1}\n".format(seq1, seq5))

print("* Testing SequenceVariation.doCopyNumberVariationSubtraction\n")

seq6 = sequenceVariationObject.doCopyNumberVariationSubtraction()

print("Copy-number variation subtraction on \n{0}\nproduces\n{1}\n".format(seq1, seq6))

print("* Testing SequenceVariation.doTransposon\n")

seq7 = sequenceVariationObject.doTransposon()

print("Transposon variation on \n{0}\nproduces\n{1}\n".format(seq1, seq7))

print("* Testing detection of mutations:\n")

variations = [seq3, seq4, seq5, seq6, seq7]
for el in variations:
    print(SequenceVariationDetection(seq1, el).detect_variation())

print("\n* Testing similarity scores:\n")

similarities = [seq2, seq3, seq4, seq5, seq6, seq7]
for el in similarities:
    print("Similarity between {0} and {1}: {2}".format(seq1, el, similarity(seq1, el)))

print("\n* Testing multiple sequence alignment:\n")

sequences = [seq1, seq2, seq3, seq4, seq5, seq6, seq7]
consensusMultipleAlign(sequences, 0.25, DNA_2)