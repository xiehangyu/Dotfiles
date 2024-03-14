from finding import *

def folded_operation(U1,U2):
    '''
    Suppose U1 and U2 are both 2-site operators, folded U2 at the bottom of U1. and the final index is 1,2,3,4 correspds to U1(1),U2(1),U1(2),U2(2), respectively.
    '''
    dim=int(np.sqrt(U1.shape[0]))
    return_matrix=np.zeros((dim**4,dim**4),dtype=complex)
    a=a+1
    for i1 in range(dim):
        for i2 in range(dim):
            for i3 in range(dim):
                for i4 in range(dim):
                    for j1 in range(dim):
                        for j2 in range(dim):
                            for j3 in range(dim):
                                for j4 in range(dim):
                                    return_matrix[i1*dim**3+i2*dim**2+i3*dim+i4,j1*dim**3+j2*dim**2+j3*dim+j4]=U1[i1*dim+i3,j1*dim+j3]*U2[i2*dim+i4,j2*dim+j4]

def depolarization_channel_folded(U,p):
    I=np.eye(2)
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    II=np.eye(4)
    XI=np.kron(X,I)
    YI=np.kron(Y,I)
    ZI=np.kron(Z,I)
    IX=np.kron(I,X)
    IY=np.kron(I,Y)
    IZ=np.kron(I,Z)
    op_list1=[(1-p)*II,p/3*XI,p/3*YI,p/3*ZI]
    op_list2=[(1-p)*II,p/3*IX,p/3*IY,p/3*IZ]
    
    
    