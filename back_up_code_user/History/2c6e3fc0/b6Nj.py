import numpy as np
from scipy.special import jv
def matrix(NA,N,t,cutoff):
    returnmatrix=np.zeros((NA,NA),dtype=complex)
    for l in range(0,NA):
        tempans=0
        for step in range(-cutoff,cutoff+1):
            tempans+=jv(step*N+l,4*t)/2
        for m1 in range(0,NA):
            if m1+l<=NA-1:
                returnmatrix[m1,m1+l]=tempans
        if l%2==0:
            tempans=tempans
        else:
            tempans=-tempans
        for m1 in range(0,NA):
            if m1-l>=0:
                returnmatrix[m1,m1-l]=tempans
    for m1 in range(0,NA):
        for m2 in range(0,NA):
            returnmatrix[m1,m2]=returnmatrix[m1,m2]*np.exp(-1j*np.pi/2*(m1+m2))
    return returnmatrix

def calculate_eigenvalue(NA,N,t,cutoff):
    returnmatrix=matrix(NA,N,t,cutoff)
    return np.linalg.eig(returnmatrix)