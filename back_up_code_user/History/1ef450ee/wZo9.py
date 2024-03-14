import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import expm


def derivative_c_channel(i,j,alpha,m,n,c_array):
    if j != n:
        return 0
    else:
        return np.conjugate(c_array[alpha,m,i])

def derivative_ccomplex_channel(i,j,alpha,m,n,c_array):
    if n != i:
        return 0
    else:
        return c_array[alpha,m,j]

def derivative_unital(i,j,alpha,m,n,c_array):
    if i!=m:
        return 0
    else:
        return np.conjugate(c_array[alpha,j,n])


def derivative_unital_complex(i,j,alpha,m,n,c_array):
    if j!=m:
        return 0
    else:
        return c_array[alpha,i,n]


def channel_condition(CM, array_c,D):
    if CM==[]:
        CM=