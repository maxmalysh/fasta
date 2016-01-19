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

