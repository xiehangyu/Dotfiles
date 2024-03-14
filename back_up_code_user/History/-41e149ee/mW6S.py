import numpy as np
from scipy.linalg import eig
from scipy.linalg import eigh
from numpy.linalg import expm

def folded_gate(gate, t):
    basis0=np.kron(np.eye(2),np.array([1,0]))
    basis1=np.kron(np.eye(2),np.array([0,1]))
    final_gate1=basis0.dot(gate.dot(basis0.T))
    final_gate2=basis0.dot(gate.dot(basis1.T))
    final_gate3=basis1.dot(gate.dot(basis0.T))
    