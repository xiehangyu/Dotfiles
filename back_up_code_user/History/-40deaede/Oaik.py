from finding import *
from CNOT_TM import dual_unitary_gate
from scipy.stats import unitary_group

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
    tildedyn4fo=tildedepchannel4fo
    key=1
    t=t-1
    while t>0:
        tildedyn4fo=np.kron(tildedyn4fo,np.eye(dim**4))
        tildedyn4fo=np.kron(np.eye(dim**4*key),tildedepchannel4fo).dot(tildedyn4fo)
        key=key+1
        t=t-1
        if t>0:
            tildedyn4fo=np.kron(tildedyn4fo,np.eye(dim**4))
            tildedyn4fo=tildedyn4fo.dot(np.kron(np.eye(dim**4*key),tildedepchannel4fo))
            key=key+1
            t=t-1
    return tildedyn4fo

    
def initial_state():
    X=np.array([[0,1],[1,0]])
    op1rank=X/np.sqrt(2)
    op1rank=np.kron(op1rank,op1rank.conj())
    return np.kron(op1rank,op1rank)

def channel_for_measurement_total(t,U,p):
    initialstate=initial_state()
    tildedyn4fo=dynamical_gate_4_folded(t,U,p)
    initialstate=np.kron(initialstate,np.eye(16**(t)))
    afterinitial=initialstate.dot(tildedyn4fo)
    I=np.eye(2)/np.sqrt(2)
    X=np.array([[0,1],[1,0]])/np.sqrt(2)
    Y=np.array([[0,-1j],[1j,0]])/np.sqrt(2)
    Z=np.array([[1,0],[0,-1]])/np.sqrt(2)
    I=np.kron(I,I)
    X=np.kron(X,X)
    Z=np.kron(Z,Z)
    Y=np.kron(Y,Y.conj())
    I=np.kron(I,I)
    X=np.kron(X,X)
    Z=np.kron(Z,Z)
    Y=np.kron(Y,Y)
    I=np.kron(np.eye(16**(t)),I)
    X=np.kron(np.eye(16**(t)),X)
    Y=np.kron(np.eye(16**(t)),Y)
    Z=np.kron(np.eye(16**(t)),Z)
    return I.dot(afterinitial)+X.dot(afterinitial)+Y.dot(afterinitial)+Z.dot(afterinitial)

def eigenvalue_for_measurement_total(t,U,p):
    ans= np.linalg.eigvals(channel_for_measurement_total(t,U,p))
    ans=sorted(ans,reverse=True)
    np.savetxt('./plot_data/eigenvalue_for_measurement_total{}and{}.txt'.format(t,p),U)
    np.savetxt('./plot_data/eigenvalue_for_measurement_total{}and{}.txt'.format(t,p),ans,append=True)
    return ans


    
def random_dual_unitary(t,p):
    intm=dual_unitary_gate(np.random.rand()*2*np.pi)
    u1=unitary_group.rvs(2)
    u2=unitary_group.rvs(2)
    eigenvalue_for_measurement_total(t,np.kron(u1,u2).dot(intm),p)