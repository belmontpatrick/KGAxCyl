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
