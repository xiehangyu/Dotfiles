import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

def period_length(i,N):
    return np.min([i,N-i])

def f_k(N,m,alpha):
    normfactor=1.0/np.power(N,1-alpha)
    k=2*np.pi/N*m
    sum_result=0
    for i in range(1,N):
        sum_result += normfactor*np.exp(-1j*k*i)/np.power(period_length(i,N)+1,alpha)
