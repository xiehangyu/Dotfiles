from finding import *

def depolarization_channel_folded(U,p):
    I=np.eye(2)
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    II=np.kron(I,I)
    XX=np.kron(X,X)
    ZZ=np.kron(Z,Z)
    YY=np.kron(Y,Y.conj())
    UU=np.kron(U,U.conj())
    
    