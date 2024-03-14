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

def distance(N,i,j):
    return min(abs(i-j),N-abs(i-j))

def construct_Hamiltonian(N,alpha,J):
    A=np.zeros([N,N])
    B=np.zeros([N,N])
    for i in range(N):
        for k in range(N):
            if i<k:
                A[i,k]=-1/(2j*np.power(N,1-alpha)*np.power(distance(N,i,k),alpha))
