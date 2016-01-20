'''
# Returns all side diagonals for any numpy matrix.
# (from upper top to the left lower corner of the matrix)
'''
import threading


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

class SynchronizedTimer:
    def __init__(self):
        self.dotplot = 0
        self.regions = 0
        self.align = 0
        self.lock = threading.Lock()

    def add(self, key, time):
        self.lock.acquire()
        current = getattr(self, key)
        setattr(self, key, current + time)
        self.lock.release()
