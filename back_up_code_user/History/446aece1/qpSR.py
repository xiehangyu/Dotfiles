import numpy as np
import matplotlib.pyplot as plt
from decomposition import kron_product
def vector_to_matrix(psi,D):
    return_matrix=np.zeros([D,D],dtype=complex)
    for i in range(D):
        for j in range(D):
            return_matrix[i,j]=psi[i*D+j]
    return return_matrix

def basis_four_identity(D):
    return kron_product([np.eye(D,dtype=complex)]*4)

def phi_plus():
    return np.array([0,1,1,0])/np.sqrt(2)

def phi_minus():
    return np.array([0,1,-1,0])/np.sqrt(2)

