import numpy as np
from scipy.special import jv
def matrix(NA,N,t,cutoff):
    returnmatrix=np.zeros(NA,NA)
    for l in range(0,NA):
        tempans=0
        for step in range(-cutoff,cutoff+1):
            tempans+=jv(step*N+l,4*t)
        for m1 in range(0,NA):
            for m2 in range(0,NA):
                returnmatrix[m1,m2]=jv((m1-m2)*N+l,4*t)/tempans 