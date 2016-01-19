from collections import namedtuple
import traceback
from pyfasta import Fasta
from Bio.SubsMat import MatrixInfo
import numpy as np
import progressbar
import matplotlib.pyplot as plt
from enum import Enum

'''
# Some helpful articles
# http://www.pnas.org/content/85/8/2444.full.pdf original one
# http://www.techfak.uni-bielefeld.de/ags/pi/lehre/KAB06/PDF/lipman_pearson_fasta.pdf
#
# http://www.gen.tcd.ie/molevol/fasta.html
# http://ab.inf.uni-tuebingen.de/teaching/ws06/albi1/script/BLAST_slides.pdf
# http://www.ch.embnet.org/CoursEMBnet/Basel03/slides/BLAST_FASTA.pdf great images
# http://vlab.amrita.edu/?sub=3&brch=274&sim=1434&cnt=1
'''

library_path = 'lib.fasta'
query_path = 'Query.txt'
library_process_limit =  None #100

smatrix = MatrixInfo.pam250
ktup = 2
best_diagonal_number = 10
band_width = 5
gap_penalty = -4

AlignmentResult = namedtuple('Result', ['title', 'score'])


'''
# Represents a diagonal run in a dot plot: a short side diagonal
# with a length not less than @ktup
'''
class Region:
    def __init__(self, x, y, length):
        self.x = x
        self.y = y
        self.length = length
        self.score = 0
        self.image = ''

    def __str__(self):
        return 'Region from [%d,%d] of length %d and score %d' % \
               (self.x, self.y, self.length, self.score)

    @property
    def x2(self):
        return self.x + self.length

    @property
    def y2(self):
        return self.y + self.length

    @property
    def diag(self):
        return self.x - self.y


class Direction(Enum):
    diag = 0
    left = 1
    up   = 2
    end  = 3


'''
# Creates a lookup table.
# It's a hash map where character corresponds to a list with it's positions in the given string
#
# "helloworld" will have:
# 'o' -> [4, 6]
# 'l' -> [2, 3, 8]
# 'h' -> [0], etc
'''
def create_lookup_table_for(seq):
    lookup_table = dict()

    for char in seq:
        lookup_table[char] = []

    for i in range(len(seq)):
        lookup_table[seq[i]].append(i)

    return lookup_table


'''
# Returns all side diagonals for any numpy matrix.
# (from upper top to the left lower corner of the matrix)
'''
def get_all_diagonals_for(matrix):
    return [
       matrix.diagonal(i) for i in range(matrix.shape[1]-1,-matrix.shape[0],-1)
    ]


def diagonal_generator(rows, cols):
    max_sum = rows + cols - 1
    for sum in range(max_sum):
        for i in range(0, rows):
            for j in range(0, cols):
                if i + j - sum == 0:
                    yield i, j

'''
# 1. Identify common words between seq1 and seq2
# 2. Score diagonals with k-word matches, identify 10 best diagonals
# 3. Rescore initial regions with a substitution score matrix
# 4. Join initial regions using gaps, penalise for gaps
# 5. Perform dynamic programming to find final alignments
'''
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

    # print(dot_plot)

    #
    # Finding diagonal regions
    # (all diagonal subsequences with length equal to or greater than ktup)
    #
    regions = []

    for k in range(-n + 1, m):
        diagonal = dot_plot.diagonal(k)
        inside = False

        for i in range(len(diagonal)):
            x, y = (i, i + k) if k >= 0 else (i - k, i)

            if not inside and diagonal[i]:
                inside = True
                x_start, y_start = x, y
                length = 1
            elif inside and not diagonal[i]:
                inside = False
                if length >= ktup:
                    regions.append(Region(x_start, y_start, length))
            elif inside:
                length += 1

        if inside and length >= ktup:
            regions.append(Region(x_start, y_start, length))

    #
    # Calculate score for each region using similarity matrix
    # Then sort by score and leave only top-10 regions
    # FIXME: find formula from the first article and "(1, 3)"
    #
    for region in regions:
        region.image = db_seq[region.x : region.x + region.length]
        region.score = sum([ smatrix[(x, x)] for x in region.image ])

    regions.sort(key=lambda run: run.score, reverse=True)
    regions = regions[:best_diagonal_number]

    #
    # Rescore those top-10 using a scoring matrix that allows
    # conservative replacements and runs of identitites shorter than ktup
    # to contribute to the similarity score
    #

    # For each of these best diagonal regions,
    # a subregion with maximal score is identified

    #
    # Join regions if this will increase score
    #
    # FIXME FIXME http://www.srmuniv.ac.in/sites/default/files/files/2(6).pdf
    #
    # FastA determines if any of the initial regions from different diagonals may be joined together to form an approximate alignment with gaps.
    # Only non-overlapping regions may be joined.
    #

    # Perform dynamic programming for survivors
    a_query, a_sequence, max_score = BoundedSmithWaterman(regions, db_seq, query_seq)
    return a_sequence, a_query, max_score


def BoundedSmithWaterman(regions, seq1, seq2):
    n, m = len(seq1), len(seq2)

    # Find bounds first
    i_min = min(regions, key=lambda reg: reg.x).x
    j_min = min(regions, key=lambda reg: reg.y).y
    i_max = max(regions, key=lambda reg: (reg.x + reg.length)).x2
    j_max = max(regions, key=lambda reg: (reg.y + reg.length)).y2

    i_min = max(i_min - band_width, 1)
    i_max = min(i_max + band_width, n)
    j_min = max(j_min - band_width, 1)
    j_max = min(j_max + band_width, m)

    left_diag = min(regions, key=lambda reg: reg.diag).diag - band_width
    right_diag = max(regions, key=lambda reg: reg.diag).diag + band_width

    # Now initialize score and traceback matrices
    max_score = 0
    scores = np.zeros(shape=(n, m), dtype=np.float)
    trace = np.empty(shape=(n, m), dtype=Direction)
    trace.fill(Direction.end)

    # Calculate scores
    for i in range(i_min, i_max):
        for j in range(j_min, j_max):
            # t - current diagonal number
            # checking whether we're inside bounds
            current_diag = i - j
            if not (current_diag >= left_diag and current_diag <= right_diag):
                continue

            pair = (seq2[j], seq1[i])
            pair = pair if pair in smatrix else (seq1[i], seq2[j])

            match  = scores[i-1][j-1] + smatrix[pair]
            insert = scores[i-1][j] + gap_penalty
            delete = scores[i][j-1] + gap_penalty
            scores[i][j] = max(match, delete, insert, 0)

            trace[(i,j)] = {
                match  : Direction.diag,
                insert : Direction.up,
                delete : Direction.left,
                0      : Direction.end,
            }[scores[i][j]]

            if scores[i][j] >= max_score:
                max_score = scores[i][j]
                i_start, j_start = i, j

    # Traceback
    aligned1 = ''
    aligned2 = ''

    i, j = i_start, j_start

    while trace[(i, j)] != Direction.end:
        direction = trace[(i, j)]

        if direction == Direction.diag:
            aligned1 = seq1[i] + aligned1
            aligned2 = seq2[j] + aligned2
            i -= 1
            j -= 1
        elif direction == Direction.up:
            aligned1 = seq1[i] + aligned1
            aligned2 = '-' + aligned2
            i -= 1
        elif direction == Direction.left:
            aligned1 = '-' + aligned1
            aligned2 = seq2[j] + aligned2
            j -= 1

    aligned1 = seq1[i] + aligned1
    aligned2 = seq2[j] + aligned2

    return aligned1, aligned2, max_score


'''
# Creates alignment outline
'''
def outline_alignment_for(aligned1, aligned2):
    alignment = ''

    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i]:
            alignment += '|'
        else:
            if aligned1[i] == '-' or aligned2[i] == '-':
                alignment += ' '
            else:
                alignment += '.'

    return alignment


'''
# Note that we should explicitly convert pyfasta record to string, as, otherwise,
# it will never be read into the memory, what will result in poor performance
# due to disk access.
'''
def process_query():
    print('Reading sequence library and query sequence')
    library = Fasta(library_path)
    queries = Fasta(query_path)
    query_sequence = str(queries["Rattus"])

    print('Processing')
    progress = progressbar.ProgressBar()
    results = []

    for record in progress(list(library.keys())[:library_process_limit]):
        try:
            library_sequence = str(library[record])
            aligned1, aligned2, score = align(library_sequence, query_sequence)
        except Exception as e:
            print(e.__class__, e)
            traceback.print_tb(e.__traceback__)
            print("Failed on "+ record)
            exit()
        results.append(AlignmentResult(title=record, score=score))

    etalone_score = sum([ smatrix[(x, x)] for x in query_sequence ])

    print("Done")
    print("Etalone score is %d" % etalone_score)
    print("Got %d results, here are top-30 among them:" % len(results))
    print("Score  | Match  | Record")

    for sequence in sorted(results, key=lambda x: x.score, reverse=True)[:30]:
        match = (sequence.score / etalone_score) * 100.0
        print("%6d | %5.3f%% | %s" % (sequence.score, match, sequence.title))

def test_simple():
    #str1 = "TACCGA"
    #str2 = "ACTGAC"

    str1 = "HEAGAWGHEE"
    str2 = "PAWHEAE"

    aligned1, aligned2, score = align(db_seq=str1, query_seq=str2)
    alignment = outline_alignment_for(aligned1, aligned2)
    print(aligned1)
    print(alignment)
    print(aligned2)
    print(score)

def test_rattus():
    library = Fasta(library_path)
    queries = Fasta(query_path)

    #query_sequence = str(queries["Rattus"])
    #query_sequence = str(library["sp|P06757|ADH1_RAT Alcohol dehydrogenase 1 OS=Rattus norvegicus GN=Adh1 PE=1 SV=3"])
    query_sequence = str(library["tr|G0Z9N0|G0Z9N0_MALDO Alcohol dehydrogenase (Fragment) OS=Malus domestica GN=Adh1-1 PE=4 SV=1"])
    library_sequence = str(library["tr|Q8K571|Q8K571_RAT ADH-like protein OS=Rattus norvegicus GN=Adh1 PE=2 SV=1"])

    aligned1, aligned2, score = align(library_sequence, query_sequence)
    outline = outline_alignment_for(aligned1, aligned2)

    with open('rattus_test.txt', 'a') as fl:
        print(query_sequence, file=fl)
        print(library_sequence, file=fl)
        print(score, file=fl)
        print(aligned1, file=fl)
        print(outline,  file=fl)
        print(aligned2, file=fl)


def main():
    process_query()
    #test_simple()
    #test_rattus()

if __name__ == "__main__":
    main()

