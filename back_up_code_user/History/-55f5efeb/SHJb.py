import numpy as np
from qutip import *
import scipy

def Gaussian_evolution(h,gamma,timeslides,N=None):
    H_hami=Qobj(h)*2
    if N is None:
        N=(np.shape(gamma)[0])//2
    rho=Qobj(gamma)/N 
    results=mesolve(H_hami,rho,timeslides,[],[])
    results_matrix=[]
    for i in range(len(timeslides)):
        results_matrix.append(results.states[i].full()*N)
    return results_matrix


def entropy_calculation(rho):
    return -np.trace(rho@scipy.linalg.logm(rho))


def construct_Hamiltonian(N,alpha,J):
    h=np.zeros((2*N,2*N))
    
