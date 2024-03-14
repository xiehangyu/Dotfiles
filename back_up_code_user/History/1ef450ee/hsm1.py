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

def derivative_unital_re(i,j,alpha,m,n,c_array):
    return derivative_unital(i,j,alpha,m,n,c_array)+derivative_unital_complex(i,j,alpha,m,n,c_array)

def derivative_unital_im(i,j,alpha,m,n,c_array):
    return 1j*(derivative_unital(i,j,alpha,m,n,c_array)-derivative_unital_complex(i,j,alpha,m,n,c_array))

def derivative_channel_re(i,j,alpha,m,n,c_array):
    return derivative_c_channel(i,j,alpha,m,n,c_array)+derivative_ccomplex_channel(i,j,alpha,m,n,c_array)

def derivative_channel_im(i,j,alpha,m,n,c_array):
    return 1j*(derivative_c_channel(i,j,alpha,m,n,c_array)-derivative_ccomplex_channel(i,j,alpha,m,n,c_array))

def channel_condition(array_c,D,d):
    CM=np.zeros([2*D**2,2*D**2*d],dtype=np.float64)
    for i in range(D):
        for j in range(D):
            for alpha in range(d):
                for m in range(D):
                    for n in range(D):
                        row_index=2*(i*D+j)
                        col_index=2*(alpha*D**2+m*D+n)
                        CM[row_index,col_index]=np.real(derivative_channel_re(i,j,alpha,m,n,array_c))
                        CM[row_index,col_index+1]=np.real(derivative_channel_im(i,j,alpha,m,n,array_c))
                        CM[row_index+1,col_index]=np.imag(derivative_channel_re(i,j,alpha,m,n,array_c))
                        CM[row_index+1,col_index+1]=np.imag(derivative_channel_im(i,j,alpha,m,n,array_c))
    return CM

def unital_condition(array_c,D,d):
    CM1=channel_condition(array_c,D,d)
    CM2=np.zeros([2*D**2,2*D**2*d],dtype=np.float64)
    for i in range(D):
        for j in range(D):
            for alpha in range(d):
                for m in range(D):
                    for n in range(D):
                        row_index=2*(i*D+j)
                        col_index=2*(alpha*D**2+m*D+n)
                        CM2[row_index,col_index]=np.real(derivative_unital_re(i,j,alpha,m,n,array_c))
                        CM2[row_index,col_index+1]=np.real(derivative_unital_im(i,j,alpha,m,n,array_c))
                        CM2[row_index+1,col_index]=np.imag(derivative_unital_re(i,j,alpha,m,n,array_c))
                        CM2[row_index+1,col_index+1]=np.imag(derivative_unital_im(i,j,alpha,m,n,array_c))
    CM_total=np.block([[CM1],[CM2]])
    return CM_total

def channel_tangent_dimension(c_array,D,d):
    CM=channel_condition(c_array,D,d)
    rank=np.linalg.matrix_rank(CM)
    print("The tangenet dimension is {}".format(2*D**2*d-rank))

def unital_tangent_dimension(c_array,D,d):
    CM=unital_condition(c_array,D,d)
    rank=np.linalg.matrix_rank(CM)
    print("The tangenet dimension is {}".format(2*D**2*d-rank))

def total_check_dimension(c_array,D,d):
    channel_tangent_dimension(c_array,D,d)
    unital_tangent_dimension(c_array,D,d)

def check_channel(c_array,D,d):
    result_matrix=np.zeros([D,D],dtype=np.complex128)
    for i in range(d):
        result_matrix+=c_array[i,:,:].conj().T@c_array[i,:,:]
    ideal=np.eye(D)
    print("The distance of channel condition is {}".format(np.linalg.norm(result_matrix-ideal)))


def unital_check(c_array,D,d):
    result_matrix=np.zeros([D,D],dtype=np.complex128)
    for i in range(d):
        result_matrix+=c_array[i,:,:]@c_array[i,:,:].conj().T
    ideal=np.eye(D)
    print("The distance of unital condition is {}".format(np.linalg.norm(result_matrix-ideal)))

def total_check(c_array,D,d):
    check_channel(c_array,D,d)
    unital_check(c_array,D,d)

def constructmethod1():
    