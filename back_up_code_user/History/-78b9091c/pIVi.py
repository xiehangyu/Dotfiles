import numpy as np
import matplotlib.pyplot as plt
from floquet_hamiltonian import foureir_transform_from_kspace_to_real_space_period1
from random_hamiltonian import calculate_entropy_from_correlation_matrix
import sys
Sampling_size=200000
def theoretical_calculation_different_order(Nsize,NA,Norder):
    N=Nsize
    if Norder==1:
        return NA**2/Nsize
    if Norder==2:
        return 2*NA**3/N**2-NA**4/N**3
    if Norder==3:
        return 6*NA**4/N**3-9*NA**5/N**4+4*NA**6/N**5
    if Norder==4:
        return 24*NA**5/N**4-72*NA**6/N**5+82*NA**7/N**6-33*NA**8/N**7
    if Norder==5:
        return 120*NA**6/N**5-600*NA**7/N**6+1250*NA**8/N**7-1225*NA**9/N**8+456*NA**10/N**9
def N_th_order_term(N,covariance_matrix):
    M=covariance_matrix
    temp=2*M-np.eye(len(M))
    temp=np.linalg.matrix_power(temp,2*N)*2
    return_value=temp.trace()
    return_value=np.real(return_value)
    return -return_value/(4*N*(2*N-1))


def generate_random_matrix_period1_in_kspace(Nsize):
    return_matrix=np.zeros([Nsize,Nsize],dtype=complex)
    for i in range(Nsize//2):
        return_matrix[i,i]=0.5
        return_matrix[i+Nsize//2,i+Nsize//2]=0.5
        phase=np.random.uniform(0,2*np.pi)
        return_matrix[i,i+Nsize//2]=0.5*np.exp(1j*phase)
        return_matrix[i+Nsize//2,i]=0.5*np.exp(-1j*phase)
    return return_matrix
    
def generate_random_matrix_period1_in_real_space(Nsize,NA):
    matrix=generate_random_matrix_period1_in_kspace(Nsize)
    UF=foureir_transform_from_kspace_to_real_space_period1(Nsize)
    real_space_matrix=UF.dot(matrix.dot(UF.T.conj()))
    return real_space_matrix[0:NA,0:NA]

def eigenvalue_entropy(x):
    if x<=0 or x>=1:
        return 0
    else:
        return -x*np.log(x)-(1-x)*np.log(1-x)

def calculate_entropy_with_CM(CM):
    eigenvalues=np.linalg.eig(CM)[0]
    eigenvalues=np.real(eigenvalues)
    sum =0
    for x in eigenvalues:
        sum=sum+eigenvalue_entropy(x)
    return sum

def Monte_Carlo_calculate(Nsize,NA,Norder):
    sum=0
    for i in range(Sampling_size):
        matrix=generate_random_matrix_period1_in_real_space(Nsize, NA)
        sum += N_th_order_term(Norder, matrix)
    return sum/Sampling_size


def different_N_fitting(NA,Norder):
    Nsize=sys.argv[1]
    Nsize=int(Nsize)
    temp=-Monte_Carlo_calculate(Nsize, NA, Norder)
    fs=open("./plot_data/trun_entropy_term_different_system_size_samplingsize{}NA{}Norder{}.txt".format(Sampling_size,NA,Norder),"a+")
    fs.write("{}\t{}\n".format(Nsize, temp))
    fs.close()

def different_NA_fitting(Nsize,Norder):
    NA=sys.argv[1]
    NA=int(NA)
    temp=-Monte_Carlo_calculate(Nsize, NA, Norder)
    fs=open("./plot_data/trun_entropy_term_different_NA_size_samplingsize{}Nsize{}Norder{}.txt".format(Sampling_size,Nsize,Norder),"a+")
    fs.write("{}\t{}\n".format(NA, temp))
    fs.close()

def different_Norder_fitting(Nsize,NA):
    Norder=sys.argv[1]
    Norder=int(Norder)
    temp=-Monte_Carlo_calculate(Nsize, NA, Norder)
    temp=temp/2*4*Norder*(2*Norder-1)
    fs=open("./plot_data/different_Norder_fittingsamplingsize{}Nsize{}NA{}.txt".format(Sampling_size,Nsize,NA),"a+")
    fs.write("{}\t{}\n".format(Norder,temp))
    fs.close()

def calculate_the_difference(Nsize,NA):
    matrix=generate_random_matrix_period1_in_real_space(Nsize, NA)
    entropy=calculate_entropy_from_correlation_matrix(matrix)
    entropy=entropy*np.log(2)
    deviations=[]
    orders=[]
    orders.append(0)
    entropy_trun=NA*np.log(2)
    deviations.append(entropy_trun-entropy)
    for i in range(1,40):
        entropy_trun_term=N_th_order_term(i, matrix)
        entropy_trun=entropy_trun+entropy_trun_term
        orders.append(i)
        deviations.append(entropy_trun-entropy)
    plt.plot(orders,deviations)
    plt.show()

if __name__=="__main__":
    different_N_fitting(80,3)