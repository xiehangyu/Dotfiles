import numpy as np
import matplotlib.pyplot as plt
from decomposition import *
import numpy.linalg as LA
def vector_to_matrix(psi,chi):
    return_matrix=np.zeros([chi,chi],dtype=complex)
    for i in range(chi):
        for j in range(chi):
            return_matrix[i,j]=psi[i*D+j]
    return return_matrix

def basis_four_identity(chi):
    return kron_product([np.eye(chi,dtype=complex)]*4)/chi**2

def phi_plus():
    return np.array([0,1,1,0])/np.sqrt(2)

def phi_minus():
    return 1j*np.array([0,-1,1,0])/np.sqrt(2)

def foldchannel_to_C_array(foldchannel,chi,D):
    whether_hermitian=LA.norm(foldchannel-foldchannel.conj().T)
    print("The norm of the difference between the foldchannel and its conjugate transpose is",whether_hermitian)
    foldchannel=(foldchannel+foldchannel.conj().T)/2
    eigenvalues,eigenvectors=LA.eigh(foldchannel)
    print(eigenvalues)
    if eigenvalues[0]<=-1e-5:
        return False
    if D<=chi**4-1 and (eigenvalues[-D-1]/eigenvalues[-1])>1e-5:
        print("The physical bond dimension is not large enough")
        return False
    else:
        C_array=np.diag(np.sqrt(eigenvalues[-D:]))@(eigenvectors[:,-D:].conj().T)
        return C_array.reshape([D,chi,chi,chi,chi])
    
def first_trial():
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    id=np.eye(2)/np.sqrt(2)
    op_x=X/2
    op_y=Y/2
    op_z=Z/2
    channel1=basis_four_identity(2)
    channel2=kron_product([id,op_x,op_x,id])    
    channel3=kron_product([op_y,op_x,id,id])
    channel4=kron_product([op_x,id,id,op_x])
    channel5=kron_product([id,id,op_y,op_x])
    return foldchannel_to_C_array(channel1+channel2+channel3+channel4+channel5,chi=2,D=16)

if __name__=='__main__':
    C_array=first_trial()
    print(decomposition_coefficient(vectorizedchannel(C_array=C_array,D=16,chi=2),chi=2))