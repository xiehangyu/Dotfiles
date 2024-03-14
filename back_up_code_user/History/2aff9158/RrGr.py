import numpy as np
from scipy.linalg import expm
def epsilonL(U,a):
    dim=np.shape(a)[0]
    temp_matrix=np.kron(np.eye(dim),a)
    temp_matrix=U.conj().transpose().dot(temp_matrix.dot(U))
    return_matrix=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim))
        return_matrix += project.dot(temp_matrix.dot(project.conj().transpose()))
    return return_matrix

def epsilonR(U,a):
    dim=np.shape(a)[0]
    temp_matrix=np.kron(a,np.eye(dim))
    temp_matrix=U.conj().transpose().dot(temp_matrix.dot(U))
    return_matrix=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim),new_vector)
        return_matrix += project.dot(temp_matrix.dot(project.conj().transpose()))
    return return_matrix

def interaction(j1,j2,j3):
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    XX=np.kron(X,X)
    YY=np.kron(Y,Y)
    ZZ=np.kron(Z,Z)
    return expm(-1j*(j1*XX+j2*YY+j3*ZZ))

def singleoperator(r,theta,phi):
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    H=r*(np.cos(theta)*Z+np.sin(theta)*cos(phi)*X+np.sin(theta)*sin(phi)*Y)
    return expm(1j*H)


def uniform_matrix(t,L,i,j,U):
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    pauli_group=[X,Y,Z]
    pa=pauli_group[i]
    pb=pauli_group[j]

