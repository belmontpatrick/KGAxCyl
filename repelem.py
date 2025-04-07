import numpy as np


def repelem(arr, repeats):
    """
    Replicates MATLAB's repelem behavior for N-dimensional arrays.
    
    Parameters:
        arr (numpy.ndarray): Input array.
        repeats (tuple): Number of repetitions along each axis.
    
    Returns:
        numpy.ndarray: Repeated array.
    """
    for axis, rep in enumerate(repeats):
        arr = np.repeat(arr, rep, axis=axis)
    return arr

# A = np.array([[1,2,3],[4,5,6],[7,8,9]])

# B = np.array([[1,2],[3,4]])

# Atile = np.tile(A, (len(B), len(B)))

# Brepelem = repelem(B, (len(A), len(A)))

# print(Atile)
# print(Brepelem)
# print('tile + repelem = ',Atile * Brepelem)

# print('kron = ', np.kron(A,B))
