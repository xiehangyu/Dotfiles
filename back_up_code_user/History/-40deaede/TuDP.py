from finding import *

def folded_operation(U1,U2)
'''
Folded U2 back to U1, suppose U1,U2 are both 2-sites operators. After the folding, the qubit 1,2,3,4 corresponds
to U1(1),U2(1),U1(2),U2(2) respectively.
'''
dim=int(np.sqrt(U1.shape[0]))
return_matrix=np.zeros((dim**4,dim**4),dtype=complex)
for(i1 in range(dim)):
    for(i2 in range(dim)):
        for(i3 in range(dim)):
            for(i4 in range(dim)):

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
    
    
    