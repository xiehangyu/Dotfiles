import numpy as np
import matplotlib.pyplot as plt
def vector_to_matrix(psi,D):
    return_matrix=np.zeros([D,D],dtype=complex)
    for i in range(D):
        for j in range(D):
            return_matrix[i,j]=psi[i*D+j]
    return return_matrix

