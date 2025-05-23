import numpy as np
from scipy.linalg import eig
from scipy.linalg import eigh
from numpy.linalg import expm

def folded_gate(gate, t):
    basis0=np.kron(np.eye(2),np.array([1,0]))
    basis1=np.kron(np.eye(2),np.array([0,1]))
    