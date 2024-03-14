import numpy as np

def R_ij_C_nmk(c_array,i,j,n,m,k):
    if j!=n:
        return 0
    else:
        return np.conjugate(c_array[i,m,k])

def R_ij_Cconj_nmk_(c_array,i,j,n,m,k):
    if i!=n:
        return 0
    else:
        return c_array[j,m,k]

def L_ij_C_nmk(c_array,i,j,n,m,k):
    if 