import numpy as np
import matplotlib.pyplot as plt
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
        print(results.states[i].tr())
        results_matrix.append(results.states[i].full()*N)
    return results_matrix


def entropy_calculation(rho):
    return -np.trace(rho@scipy.linalg.logm(rho))

def distance(N,i,j):
    return min(abs(i-j),N-abs(i-j))

def construct_Hamiltonian(N,alpha,J):
    A=np.zeros([N,N],dtype=complex)
    B=np.zeros([N,N],dtype=complex)
    for i in range(N):
        for k in range(N):
            if i<k:
                A[i,k]=-1/(2j*np.power(N,1-alpha)*np.power(distance(N,i,k),alpha))
                B[i,k]=-1/(2j*np.power(N,1-alpha)*np.power(distance(N,i,k),alpha))
            elif i==k:
                A[i,k]=-J/2.0
                B[i,k]=0
            elif i>k:
                A[i,k]=1/(2j*np.power(N,1-alpha)*np.power(distance(N,i,k),alpha))
                B[i,k]=1/(2j*np.power(N,1-alpha)*np.power(distance(N,i,k),alpha))
    h=np.block([[-np.conj(A),B],[-np.conj(B),A]])
    return h

def simple_initial_state(N):
    one_particle=np.array([[1,0],[0,0]])
    return np.kron(one_particle,np.eye(N))


def simple_initial_state_entropy_evolution(N,alpha,J,total_time,dt=0.2):
    gamma=simple_initial_state(N)
    h=construct_Hamiltonian(N,alpha,J)
    timeslides=np.arange(0,total_time,dt)
    result_gammas=Gaussian_evolution(h,gamma,timeslides,N)
    entropys=[]
    for i in range(len(result_gammas)):
        entropys.append(entropy_calculation(result_gammas[i]))
    plt.plot(timeslides,entropys)
    plt.xlabel('time')
    plt.ylabel('entropy')
    plt.show()
    return entropys
