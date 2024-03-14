import numpy as np
import matplotlib.pyplot as plt
from decomposition import kron_product
import numpy.linalg as LA
def vector_to_matrix(psi,chi):
    return_matrix=np.zeros([chi,chi],dtype=complex)
    for i in range(chi):
        for j in range(chi):
            return_matrix[i,j]=psi[i*D+j]
    return return_matrix

def basis_four_identity(chi):
    return kron_product([np.eye(chi,dtype=complex)]*4)

def phi_plus():
    return np.array([0,1,1,0])/np.sqrt(2)

def phi_minus():
    return 1j*np.array([0,-1,1,0])/np.sqrt(2)

def foldchannel_to_C_array(foldchannel,chi,D):
    foldchannel=foldchannel.reshape([chi]*8).transpose([0,2,4,6,1,3,5,7]).reshape([chi**4]*2)
    whether_hermitian=LA.norm(foldchannel-foldchannel.conj().T)
    print("The norm of the difference between the foldchannel and its conjugate transpose is",whether_hermitian)
    
    