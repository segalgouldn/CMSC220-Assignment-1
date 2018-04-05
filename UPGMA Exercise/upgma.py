# Noah Segal-Gould, Aaron Krapf, and Zoe Terhune
# 4/4/18
# The three of us collaborated on this assignment.
# We made use of the textbook as well as code on Moodle 2.
# For the example in the original code that was on Moodle 2,
# The resultant tree is different now that we've implemented UPGMA.

import numpy as np

STANDARD_GENETIC_CODE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UCU': 'Ser', 'UCC': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UGU': 'Cys', 'UGC': 'Cys',
    'UUA': 'Leu', 'UCA': 'Ser', 'UAA': None, 'UGA': None,
    'UUG': 'Leu', 'UCG': 'Ser', 'UAG': None, 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CCU': 'Pro', 'CCC': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CGU': 'Arg', 'CGC': 'Arg',
    'CUA': 'Leu', 'CUG': 'Leu', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAA': 'Gln', 'CAG': 'Gln', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'ACU': 'Thr', 'ACC': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AGU': 'Ser', 'AGC': 'Ser',
    'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
    'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GCU': 'Ala', 'GCC': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GGU': 'Gly', 'GGC': 'Gly',
    'GUA': 'Val', 'GUG': 'Val', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAA': 'Glu', 'GAG': 'Glu', 'GGA': 'Gly', 'GGG': 'Gly'}

GES_SCALE = {'F': -3.7, 'M': -3.4, 'I': -3.1, 'L': -2.8, 'V': -2.6,
             'C': -2.0, 'W': -1.9, 'A': -1.6, 'T': -1.2, 'G': -1.0,
             'S': -0.6, 'P': 0.2, 'Y': 0.7, 'H': 3.0, 'Q': 4.1,
             'N': 4.8, 'E': 8.2, 'K': 8.8, 'D': 9.2, 'R': 12.3}

DNA_1 = {'G': {'G': 1, 'C': 0, 'A': 0, 'T': 0},
         'C': {'G': 0, 'C': 1, 'A': 0, 'T': 0},
         'A': {'G': 0, 'C': 0, 'A': 1, 'T': 0},
         'T': {'G': 0, 'C': 0, 'A': 0, 'T': 1}}

REV_COMP = {'G': {'G': -1, 'C': 1, 'A': -1, 'T': -1},
            'C': {'G': 1, 'C': -1, 'A': -1, 'T': -1},
            'A': {'G': -1, 'C': -1, 'A': -1, 'T': 1},
            'T': {'G': -1, 'C': -1, 'A': 1, 'T': -1}}

DNA_2 = {'G': {'G': 1, 'C': -3, 'A': -3, 'T': -3, 'N': 0},
         'C': {'G': -3, 'C': 1, 'A': -3, 'T': -3, 'N': 0},
         'A': {'G': -3, 'C': -3, 'A': 1, 'T': -3, 'N': 0},
         'T': {'G': -3, 'C': -3, 'A': -3, 'T': 1, 'N': 0},
         'N': {'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N': 0}}

BLOSUM62 = {'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1,
                  'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'X': 0},
            'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3,
                  'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'X': 0},
            'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3,
                  'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'X': 0},
            'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3,
                  'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'X': 0},
            'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1,
                  'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'X': 0},
            'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3,
                  'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'X': 0},
            'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3,
                  'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'X': 0},
            'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4,
                  'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'X': 0},
            'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3,
                  'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'X': 0},
            'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4,
                  'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'X': 0},
            'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2,
                  'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'X': 0},
            'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3,
                  'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'X': 0},
            'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1,
                  'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'X': 0},
            'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0,
                  'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'X': 0},
            'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3,
                  'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'X': 0},
            'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2,
                  'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'X': 0},
            'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1,
                  'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'X': 0},
            'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3,
                  'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'X': 0},
            'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1,
                  'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'X': 0},
            'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3,
                  'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'X': 0},
            'X': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
                  'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0, 'X': 0}}


def profile(alignment):
    n = len(alignment[0])
    nSeq = float(len(alignment))
    prof = []

    for i in range(n):
        counts = {}

        for seq in alignment:
            letter = seq[i]
            if letter == '-':
                continue

            counts[letter] = counts.get(letter, 0) + 1

        for letter in counts:
            counts[letter] /= nSeq

        prof.append(counts)

    return prof


def profileAlign(profileA, profileB, simMatrix, insert=8, extend=4):
    numI = len(profileA) + 1
    numJ = len(profileB) + 1

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

            if routeMatrix[i - 1][j] == 1:
                penalty1 = extend

            elif routeMatrix[i][j - 1] == 2:
                penalty2 = extend

            fractionsA = profileA[i - 1]
            fractionsB = profileB[j - 1]

            similarity = 0.0
            totalWeight = 0.0
            for residueA in fractionsA:
                for residueB in fractionsB:
                    weight = fractionsA[residueA] * fractionsB[residueB]
                    totalWeight += weight
                    similarity += weight * simMatrix[residueA][residueB]

            penalty1 *= totalWeight
            penalty2 *= totalWeight

            paths = [scoreMatrix[i - 1][j - 1] + similarity,  # Route 0
                     scoreMatrix[i - 1][j] - penalty1,  # Route 1
                     scoreMatrix[i][j - 1] - penalty2]  # Route 2

            best = max(paths)
            route = paths.index(best)

            scoreMatrix[i][j] = best
            routeMatrix[i][j] = route

    profileOutA = []
    profileOutB = []

    i = numI - 1
    j = numJ - 1
    score = scoreMatrix[i][j]

    while i > 0 or j > 0:
        route = routeMatrix[i][j]

        if route == 0:  # Diagonal
            profileOutA.append(profileA[i - 1])
            profileOutB.append(profileB[j - 1])
            i -= 1
            j -= 1

        elif route == 1:  # Gap in profile B
            profileOutA.append(profileA[i - 1])
            profileOutB.append(None)
            i -= 1

        elif route == 2:  # Gap in profile A
            profileOutA.append(None)
            profileOutB.append(profileB[j - 1])
            j -= 1

    profileOutA.reverse()
    profileOutB.reverse()

    return score, profileOutA, profileOutB


def calcSeqSimilarity(seqA, seqB, simMatrix):
    numPlaces = min(len(seqA), len(seqB))

    totalScore = 0.0

    for i in range(numPlaces):
        residueA = seqA[i]
        residueB = seqB[i]

        totalScore += simMatrix[residueA][residueB]

    return totalScore


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

            if routeMatrix[i - 1][j] == 1:
                penalty1 = extend

            elif routeMatrix[i][j - 1] == 2:
                penalty2 = extend

            similarity = simMatrix[seqA[i - 1]][seqB[j - 1]]

            paths = [scoreMatrix[i - 1][j - 1] + similarity,  # Route 0
                     scoreMatrix[i - 1][j] - penalty1,  # Route 1
                     scoreMatrix[i][j - 1] - penalty2]  # Route 2

            best = max(paths)
            route = paths.index(best)

            scoreMatrix[i][j] = best
            routeMatrix[i][j] = route

    alignA = []
    alignB = []

    i = numI - 1
    j = numJ - 1
    score = scoreMatrix[i][j]

    while i > 0 or j > 0:
        route = routeMatrix[i][j]

        if route == 0:  # Diagonal
            alignA.append(seqA[i - 1])
            alignB.append(seqB[j - 1])
            i -= 1
            j -= 1

        elif route == 1:  # Gap in seqB
            alignA.append(seqA[i - 1])
            alignB.append('-')
            i -= 1

        elif route == 2:  # Gap in seqA
            alignA.append('-')
            alignB.append(seqB[j - 1])
            j -= 1

    alignA.reverse()
    alignB.reverse()
    alignA = ''.join(alignA)
    alignB = ''.join(alignB)

    return score, alignA, alignB


def getDistanceMatrix(seqs, simMatrix):
    n = len(seqs)
    # NbyN
    matrix = [[0.0] * n for _ in range(n)]

    maxScores = [calcSeqSimilarity(x, x, simMatrix) for x in seqs]

    for i in range(n - 1):
        seqA = seqs[i]

        for j in range(i + 1, n):
            seqB = seqs[j]

            score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix)
            maxScore = max(maxScores[i], maxScores[j])
            dist = maxScore - score

            matrix[i][j] = dist
            matrix[j][i] = dist

    return matrix


def getDistToJunction(distMatrix, i, j):
    n = len(distMatrix)
    row = distMatrix[i]
    column = distMatrix[j]

    dist = distMatrix[i][j] + (sum(row) - sum(column)) / (n - 2)
    dist *= 0.5

    return dist


def getJoinPair(distMatrix):
    n = len(distMatrix)
    minQ = None
    joinPair = None

    for i in range(n - 1):
        #
        sumRow = sum(distMatrix[i])

        #
        for j in range(i + 1, n):
            #
            sumNext = sum(distMatrix[j])

            #
            dist = distMatrix[i][j]
            q = (n - 2) * dist - sumRow - sumNext
            print(i, j, q)

            if (minQ is None) or (q < minQ):
                minQ = q
                joinPair = [i, j]

    return joinPair


def neighbourJoinTree(distMatrix):
    joinOrder = []
    n = len(distMatrix)
    tree = list(range(n))

    while n > 2:

        x, y = getJoinPair(distMatrix)

        node = (tree[x], tree[y])
        joinOrder.append(node)
        tree.append(node)

        del tree[y]
        del tree[x]

        distX = getDistToJunction(distMatrix, x, y)
        distY = getDistToJunction(distMatrix, y, x)

        distMatrix.append([0] * (n + 1))

        for i in range(n):
            if i not in (x, y):
                dist = (distMatrix[x][i] - distX) + (distMatrix[y][i] - distY)
                dist *= 0.5

                distMatrix[i].append(dist)
                distMatrix[n][i] = dist

        del distMatrix[y]
        del distMatrix[x]

        for row in distMatrix:
            del row[y]
            del row[x]

        n -= 1

    tree = tuple(tree)
    joinOrder.append(tree)

    return tree, joinOrder


def show_result(result):
    for i in result:
        print("-> Combined: ({0}, {1})".format(i[0], i[1]))
    print()


def get_min_value(matrix):
    min_val = float("inf")
    pos = (-1, -1)
    length = len(matrix)
    for i in range(length):
        for j in range(length):
            if matrix[i][j] < min_val and matrix[i][j] != 0:
                min_val = matrix[i][j]
                pos = (i, j)
    return min_val, pos


def show_matrix(labels, matrix):
    row_length = len(matrix[0])
    row_string = " " * 25 + " ".join("{:>24}" for _ in range(row_length)).format(*labels)
    print("{0}\n{1}".format(row_string, "-" * len(row_string)))
    for label, row in zip(labels, matrix):
        print("{:>25}{}".format(label, " ".join("{:>24}" for _ in range(row_length)).format(*row)))
    print("-" * len(row_string) + "\n")

def run_tests():
    seqs = ['QPVHPFSRPAPVVIILIILCVMAGVIGTILLISYGIRLLIK',
            'QLVHRFTVPAPVVIILIILCVMAGIIGTILLISYTIRRLIK',
            'QLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKKSPSDVKPLPSPD',
            'QLVHEFSELVIALIIFGVMAGVIGTILFISYGSRRLIKKSESDVQPLPPPD',
            'MLEHEFSAPVAILIILGVMAGIIGIILLISYSIGQIIKKRSVDIQPPEDED',
            'PIQHDFPALVMILIILGVMAGIIGTILLISYCISRMTKKSSVDIQSPEGGD',
            'QLVHIFSEPVIIGIIYAVMLGIIITILSIAFCIGQLTKKSSLPAQVASPED',
            'LAHDFSQPVITVIILGVMAGIIGIILLLAYVSRRLRKRPPADVP',
            'SYHQDFSHAEITGIIFAVMAGLLLIIFLIAYLIRRMIKKPLPVPKPQDSPD']

    matrix = getDistanceMatrix(seqs, BLOSUM62)

    result = []
    length = len(matrix)
    labels = [str(i) for i in range(length)]
    print("Start:")
    show_matrix(labels, matrix)
    step = 1
    while length > 2:
        value, pos = get_min_value(matrix)
        for i in range(length):
            if i != pos[0]:
                matrix[i][pos[0]] = round((matrix[i][pos[0]] + matrix[i][pos[1]]) / 2.0, 2)
                matrix[pos[0]][i] = matrix[i][pos[0]]
        matrix.pop(pos[1])
        for row in matrix:
            row.pop(pos[1])
        result.append((labels[pos[0]], labels[pos[1]]))
        insert = labels[pos[0]], labels[pos[1]]
        if len(insert[0]) == 1 and len(insert[1]) == 1:
            labels[pos[0]] = "{0}{1}".format(insert[0], insert[1])
        else:
            labels[pos[0]] = "({0}, {1})".format(insert[0], insert[1])
        labels.pop(pos[1])
        print("Step: {}".format(step))
        show_matrix(labels, matrix)
        step += 1
        show_result(result)
        result = []
        length -= 1
    print("Result should be: {}".format("(((7, (01)), (45)), ((23), (68))"))
    
if __name__ == "__main__":
    run_tests()
    