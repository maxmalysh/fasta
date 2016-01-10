# 1. Identify common words between seq1 and seq2
# 2. Score diagonals with k-word matches, identify 10 best diagonals
# 3. Rescore initial regions with a substitution score matrix
# 4. Join initial regions using gaps, penalise for gaps
# 5. Perform dynamic programming to find final alignments
from collections import namedtuple
from pyfasta import Fasta
from Bio.SubsMat import MatrixInfo
import numpy as np
import progressbar
import matplotlib.pyplot as plt


library_path = 'lib.fasta'
query_path = 'Query.txt'

smatrix = MatrixInfo.pam250 # MatrixInfo.blosum50
ktup = 2
best_diagonals = 10
gap_open = -10
gap_extend = -2

#
# Nice articles
# http://www.gen.tcd.ie/molevol/fasta.html
# http://ab.inf.uni-tuebingen.de/teaching/ws06/albi1/script/BLAST_slides.pdf
# http://www.ch.embnet.org/CoursEMBnet/Basel03/slides/BLAST_FASTA.pdf
# http://vlab.amrita.edu/?sub=3&brch=274&sim=1434&cnt=1
#

AlignmentResult = namedtuple('Result', ['sequence', 'score'])


#
# Represents a diagonal run in a dot plot: a short side diagonal
# with a length not less than @ktup
#
class DiagonalRun:
    def __init__(self, x, y, length):
        self.x = x
        self.y = y
        self.length = length
        self.score = 0

    def score_from(self, sequence):
        ...

#
# Creates a lookup table.
# It's a hash map where character corresponds to a list with it's positions in the give string
#
# "helloworld" will have:
# 'o' -> [4, 6]
# 'l' -> [2, 3, 8]
# 'h' -> [0], etc
#
def create_lookup_table_for(seq):
    lookup_table = dict()

    for char in seq:
        lookup_table[char] = []

    for i in range(len(seq)):
        lookup_table[seq[i]].append(i)

    return lookup_table


#
# Returns all side diagonals for any numpy matrix.
# (from upper top to the left lower corner of the matrix)
#
def get_all_diagonals_for(matrix):
    return [
       matrix.diagonal(i) for i in range(matrix.shape[1]-1,-matrix.shape[0],-1)
    ]

def align(db_seq, query_seq):
    n, m = len(db_seq), len(query_seq)

    # Create a lookup table
    seq2_lookup_table = create_lookup_table_for(query_seq)

    # Build a dot plot
    dot_plot = np.zeros(shape=(n, m), dtype=np.int)
    for i in range(n):
        if db_seq[i] not in query_seq:
            continue
        for j in range(m):
            if j in seq2_lookup_table[db_seq[i]]:
                dot_plot[i,j] = 1

    print(dot_plot)

    #
    # Finding hot-spots
    # (all diagonal subsequences with length equal to or greater than ktup)
    #

    diagonals = get_all_diagonals_for(dot_plot)
    ktup_runs = []

    for k in range(len(diagonals)):
        diagonal = diagonals[k]
        inside = False

        for i in range(len(diagonal)):
            if not inside and diagonal[i]:
                inside = True
                x_start = n - (len(diagonal) - 1) + i
                y_start = i
                length = 1
            elif inside and not diagonal[i]:
                inside = False
                if length >= ktup:
                    ktup_runs.append(DiagonalRun(x_start, y_start, length))
            elif inside:
                length += 1
        if inside:
            ktup_runs.append(DiagonalRun(x_start, y_start, length))

        print(ktup_runs)



    # Apply smith-waterman



    aligned1, aligned2 = db_seq, query_seq
    score = 0
    return aligned1, aligned2, score




def outline_alignment_for(seq1, seq2):
    outline = ''

    for char1, char2 in zip(seq1, seq2):
        if char1 == char2:
            outline += char1
        elif char1 == '-' or char2 == '-':
            outline += ' '
        else:
            outline += '·'

    return outline

def process_query():
    print('Reading sequence library and query sequence')
    library = Fasta(library_path)
    queries = Fasta(query_path)
    query_sequence = queries["Rattus"]

    print('Processing')
    progress = progressbar.ProgressBar()
    results = []

    for record in progress(library.keys()):
        library_sequence = library[record]
        aligned1, aligned2, score = align(library_sequence, query_sequence)
        results.append(AlignmentResult(sequence=library_sequence, score=score))

def main():

    str1 = "TACCGA"
    str2 = "ACTGAC"

    aligned1, aligned2, _ = align(str1, str2)
    #process_query()


    # Match input with each sequence from the library
    # For each sequence:
    # 1. We're aligning input sequence and a library sequence, getting a similarity score
    # 2.



        # time.sleep(2.0/len(library))
        # if len(library_sequence) <= 48:
        #     results[library_sequence] = len(library_sequence) / 48.0 * 100

    #
    # Print results
    #
    # print("Done")
    # print("Got %d results, here are top-10 among them:" % len(results))
    # print("Score  | Record")
    # for sequence in sorted(results, key=results.get, reverse=True):
    #     print("%6.2f | %s" % (results[sequence], sequence))

if __name__ == "__main__":
    main()

