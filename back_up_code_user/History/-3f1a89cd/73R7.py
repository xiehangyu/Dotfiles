import numpy as np
import matplotlib.pyplot as plt

def operatornorm(A):
    """Compute the operator norm of a matrix A."""
    return np.linalg.norm(A, ord=2)

def distance_PC(i,j,N):
    """Compute the distance between two sites on a periodic chain."""
    return min(abs(i-j),N-abs(i-j))


def norm_fourier_op(N,m):
    """Compute the norm of the Fourier operator."""
    k=2*np.pi*m/N
    Z=np.array([[1,0],[0,-1]])
    O=np.zeros([2**N,2**N])
    for i in range(N):
        O=O+np.exp(1j*k*i)*np.kron(np.kron(np.eye(2**(i)),Z),np.eye(2**(N-i-1)))
    return operatornorm(O)/N


def norm_fourier_op_with_frequency(N):
    ms=[]
    norms=[]
    for m in range(N):
        ms.append(m)
        norms.append(norm_fourier_op(N,m))
    plt.plot(ms,norms)
    plt.show()
            
def norm_fourier_op_with_N(m=1):
    Ns=[]
    norms=[]
    for N in range(3,13):
        Ns.append(N)
        norms.append(norm_fourier_op(N,m))
    plt.plot(Ns,norms)
    plt.show()


def Fourier_longrange_transformation(N,m,alpha):
    k=2*np.pi*m/N
    returnvalue=0
    for i in range(1,N):
        returnvalue += np.exp(-1j*k*i)/distance_PC(i,0,N)**alpha
    return returnvalue


def Fourier_longrange_with_m(N,alpha):
    ms=[]
    norms=[]
    for m in range(1,N):
        ms.append(np.abs(m))
        norms.append(Fourier_longrange_transformation(N,m,alpha))
    plt.plot(ms,norms)
    plt.show()