from finding import *

def optimized_function(U1, U2):
    dim = int(np.sqrt(U1.shape[0]))
    U1_reshaped = U1.reshape(dim, dim, dim, dim)
    U2_reshaped = U2.reshape(dim, dim, dim, dim)
    return_matrix = np.einsum('acxy,bdwz->abcdxwyz', U1_reshaped, U2_reshaped).reshape(dim**4, dim**4)
    return return_matrix

def folded_operation(U1,U2):
    '''
    Suppose U1 and U2 are both 2-site operators, folded U2 at the bottom of U1. and the final index is 1,2,3,4 correspds to U1(1),U2(1),U1(2),U2(2), respectively.
    '''
    dim=int(np.sqrt(U1.shape[0]))
    return_matrix=np.zeros((dim**4,dim**4),dtype=complex)
    for i1 in range(dim):
        for i2 in range(dim):
            for i3 in range(dim):
                for i4 in range(dim):
                    for j1 in range(dim):
                        for j2 in range(dim):
                            for j3 in range(dim):
                                for j4 in range(dim):
                                    return_matrix[i1*dim**3+i2*dim**2+i3*dim+i4,j1*dim**3+j2*dim**2+j3*dim+j4]=U1[i1*dim+i3,j1*dim+j3]*U2[i2*dim+i4,j2*dim+j4]
    return return_matrix



def depolarization_channel_folded(U,p):
    I=np.eye(2)
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    II=np.eye(4)
    XX=np.kron(X,X)
    ZZ=np.kron(Z,Z)
    YY=np.kron(Y,Y.conj())
    op_list1=[(1-p)*np.eye(16),p/3*np.kron(XX,II),p/3*np.kron(YY,II),p/3*np.kron(ZZ,II)]
    op_list2=[(1-p)*np.eye(16),p/3*np.kron(II,XX),p/3*np.kron(II,YY),p/3*np.kron(II,ZZ)]
    U_comp=U.conj()
    temp_matrix=folded_operation(U,U_comp)
    return_matrix=np.zeros(np.shape(temp_matrix),dtype=complex)
    for i in op_list1:
        for j in op_list2:
            return_matrix+=i.dot(j).dot(temp_matrix)
    return return_matrix


def dynamical_gate_4_folded(t,U,p):
    depchannel2fo=depolarization_channel_folded(U,p)
    tildedepchannel2fo=generaltildeU(depchannel2fo)
    tildedepchannel4fo=optimized_function(tildedepchannel2fo,tildedepchannel2fo)
    dim=int(np.sqrt(U.shape[0]))
        

    
    
    