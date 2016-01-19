import numpy as np
import matplotlib.pyplot as plt
from Bio.SubsMat import MatrixInfo
from enum import Enum
import time

'''
# Some helpful articles
# http://www.pnas.org/content/85/8/2444.full.pdf original one
# http://www.techfak.uni-bielefeld.de/ags/pi/lehre/KAB06/PDF/lipman_pearson_fasta.pdf
#
# http://www.gen.tcd.ie/molevol/fasta.html
# http://ab.inf.uni-tuebingen.de/teaching/ws06/albi1/script/BLAST_slides.pdf
# http://www.ch.embnet.org/CoursEMBnet/Basel03/slides/BLAST_FASTA.pdf great images
# http://vlab.amrita.edu/?sub=3&brch=274&sim=1434&cnt=1
# http://www.cs.tau.ac.il/~rshamir/algmb/01/scribe03/lec03.pdf – best description on how to MERGE regions
#
'''

smatrix = MatrixInfo.pam250
ktup = 2
best_diagonal_number = 10
band_width = 8
gap_penalty = -4


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


'''
# Creates a lookup table.
# It's a hash map where character corresponds to a list with it's positions in the given string
#
# "helloworld" will have:
# 'o' -> [4, 6]
# 'l' -> [2, 3, 8]
# 'h' -> [0], etc
'''

def get_lookup_table_for(seq):
    lookup_table = {}

    for char in seq:
        lookup_table[char] = []

    for i in range(len(seq)):
        lookup_table[seq[i]].append(i)

    return lookup_table

dotplot_time = 0
region_time = 0
align_time = 0

def get_performance():
    return dotplot_time, region_time, align_time

'''
# 1. Identify common words between seq1 and seq2
# 2. Score diagonals with k-word matches, identify 10 best diagonals
# 3. Rescore initial regions with a substitution score matrix
# 4. Join initial regions using gaps, penalise for gaps
# 5. Perform dynamic programming to find final alignments
'''
def align(db_seq, query_seq):
    global dotplot_time, region_time, align_time

    n, m = len(db_seq), len(query_seq)

    # Create a lookup table
    db_lookup_table = get_lookup_table_for(db_seq)

    # Build a dot plot
    t1 = time.perf_counter()

    dot_plot = np.zeros(shape=(n,m), dtype=np.int)
    for j in range(m):
        char = query_seq[j]
        if char in db_lookup_table:
            for i in db_lookup_table[char]:
                dot_plot[i, j] = 1

    dotplot_time += time.perf_counter() - t1

    #
    # Finding diagonal regions
    # (all diagonal subsequences with length equal to or greater than ktup)
    #
    t1 = time.perf_counter()
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

    region_time += time.perf_counter() - t1

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
    t1 = time.perf_counter()
    if len(regions) != 0:
        best_region = regions[0]
        db_aligned, query_aligned, max_score = BoundedSmithWaterman(regions[:3], db_seq, query_seq)
    else:
        db_aligned, query_aligned, max_score = '', '', 0
    align_time += time.perf_counter() - t1

    return db_aligned, query_aligned, max_score


class Direction(Enum):
    diag = 0
    left = 1
    up   = 2
    end  = 3

def BoundedSmithWaterman(regions, seq1, seq2):
    n, m = len(seq1), len(seq2)

    # Find bounds first
    i_min = min(regions, key=lambda reg: reg.x).x
    j_min = min(regions, key=lambda reg: reg.y).y
    i_max = max(regions, key=lambda reg: reg.x2).x2
    j_max = max(regions, key=lambda reg: reg.y2).y2

    i_min = max(i_min - band_width, 1)
    i_max = min(i_max + band_width, n)
    j_min = max(j_min - band_width, 1)
    j_max = min(j_max + band_width, m)

    left_diag = min(regions, key=lambda reg: reg.diag).diag - band_width
    right_diag = max(regions, key=lambda reg: reg.diag).diag + band_width


    # i_min = max(1, region.x  - band_width)
    # i_max = min(n, region.x2 + band_width)
    # j_min = max(1, region.y  - band_width)
    # j_max = min(m, region.y2 + band_width)
    #
    # left_diag  = region.diag - band_width
    # right_diag = region.diag + band_width

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
            if not (current_diag > left_diag and current_diag < right_diag):
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

