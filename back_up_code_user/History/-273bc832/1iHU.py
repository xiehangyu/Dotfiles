import numpy as np
import scipy
import cvxpy as cp
from copy import deepcopy
from ncon import ncon
from decomposition import kron_product

def Gn(Oas,varphis,c_a,P):
    number_oas=len(Oas)
    number_phis=len(varphis)
    return_array=np.zeros([number_oas+number_phis**2],dtype=complex)
    for i in range(number_oas):
        temp_value=0
        for j in range(number_oas):
            temp_value += np.trace(Oas[i].dot(Oas[j]))*c_a[j]
        for alpha in range(number_phis):
            for beta in range(number_phis):
                temp_value -=np.conj(varphis[beta]).T.dot(Oas[i]).dot(varphis[alpha])*P[alpha,beta]
                #temp_value +=np.conj(varphis[alpha]).T.dot(Oas[i]).dot(varphis[beta])*P[alpha,beta]
        return_array[i]=temp_value
    for gamma in range(number_phis):
        for delta in range(number_phis):
            temp_value=0
            for i in range(number_oas):
                temp_value -= np.conj(varphis[gamma]).T.dot(Oas[i]).dot(varphis[delta])*c_a[i]
                #temp_value += np.conj(varphis[delta]).T.dot(Oas[i]).dot(varphis[gamma])*c_a[i]
            temp_value += P[gamma,delta]
            return_array[number_oas+gamma*number_phis+delta]=temp_value
    return return_array

def construct_M_R_matrix(Oas,varphis):
    M=np.zeros([len(Oas),len(Oas)],dtype=complex)
    R=np.zeros([len(Oas),len(varphis)**2],dtype=complex)
    for i in range(len(Oas)):
        for j in range(len(Oas)):
            M[i,j]=np.trace(Oas[i].dot(Oas[j]))
    state_matrix=np.hstack(varphis)
    for i in range(len(Oas)):
        print(i,flush=True)
        temp_matrix=-np.conj(state_matrix).T.dot(Oas[i]).dot(state_matrix)
        for alpha in range(len(varphis)):
            for beta in range(len(varphis)):
                R[i,alpha*len(varphis)+beta]=temp_matrix[beta,alpha]
            #R[i,alpha*len(varphis)+beta]=-np.conj(varphis[beta]).T.dot(Oas[i]).dot(varphis[alpha])
    return M,R

def matrix_multiplication(M,R,X):
    Y=np.zeros(len(X),dtype=complex)
    Y[:len(M)]=np.block([M,R]).dot(X)
    Y[len(M):]+=X[len(M):]
    Y[len(M):]+=np.conj(R).T.dot(X[:len(M)])
    return Y


def Optimization_by_hand_V2(Oas,varphis):
    M,R=construct_M_R_matrix(Oas,varphis)
    c_a=np.zeros(len(Oas),dtype=complex)
    P=3*np.eye(len(varphis),dtype=complex)
    X=np.concatenate((c_a,P.flatten()))
    Xlast=np.zeros(len(X),dtype=complex)
    while np.linalg.norm(X-Xlast)>5e-7:
        G=matrix_multiplication(M,R,X)
        eta=np.vdot(G,G)/np.vdot(G,matrix_multiplication(M,R,G))
        print(eta,flush=True)
        Xlast=X
        X=X-eta*G
        c_a=X[:len(Oas)]
        P=X[len(Oas):].reshape([len(varphis),len(varphis)])
        T,U=scipy.linalg.schur(P)
        np.fill_diagonal(T, np.where(np.diag(T) <= 0.1, 0.1, np.diag(T)))
        P=U@T@U.conj().T
        first_term=sum(c_a[i]*Oas[i] for i in range(len(Oas)))
        second_term=sum(P[i,j]*np.kron(varphis[i],np.conj(varphis[j]).T) for i in range(len(varphis)) for j in range(len(varphis)))
    first_term=sum(c_a[i]*Oas[i] for i in range(len(Oas)))
    second_term=sum(P[i,j]*np.kron(varphis[i],np.conj(varphis[j]).T) for i in range(len(varphis)) for j in range(len(varphis)))
    lambda_=np.linalg.norm(first_term-second_term)
    return c_a,P,lambda_


def test_found_operator(c_a,Oas,varphis,M1):
    first_term=sum(c_a[i]*Oas[i] for i in range(len(Oas)))
    print(np.linalg.norm(first_term@M1),flush=True)
    orthogonal_matrix=np.hstack(varphis)
    projection_to_orthogonal=np.conj(orthogonal_matrix).T.dot(first_term).dot(orthogonal_matrix)
    print(scipy.linalg.eigvalsh(projection_to_orthogonal),flush=True)


def Optimization_by_hand(Oas,varphis):
    c_a=np.zeros(len(Oas),dtype=complex)
    P=np.eye(len(varphis),dtype=complex)
    X=np.concatenate((c_a,P.flatten()))
    eta=1
    Xlast=np.zeros(len(X),dtype=complex)
    Glast=np.zeros(len(X),dtype=complex)
    while np.linalg.norm(X-Xlast)>5e-7:
        G=Gn(Oas,varphis,c_a,P)
        deltaX=X-Xlast
        deltaG=G-Glast
        eta=np.vdot(deltaX,deltaX)/np.vdot(deltaX,deltaG)
        print("eta=",eta)
        #eta=np.min([eta,1/len(varphis)])
        #if eta==-np.inf or eta==np.inf or np.isnan(eta):
        #    eta=1
        Xlast=X
        Glast=G
        X=X-2*eta*G
        c_a=X[:len(Oas)]
        print(c_a)
        P=X[len(Oas):].reshape([len(varphis),len(varphis)])
        T,U=scipy.linalg.schur(P)
        np.fill_diagonal(T, np.where(np.diag(T) <= 1, 1, np.diag(T)))
        P=U@T@U.conj().T
        first_term=sum(c_a[i]*Oas[i] for i in range(len(Oas)))
        second_term=sum(P[i,j]*np.kron(varphis[i],np.conj(varphis[j]).T) for i in range(len(varphis)) for j in range(len(varphis)))
        print(np.linalg.norm(first_term-second_term))
    first_term=sum(c_a[i]*Oas[i] for i in range(len(Oas)))
    second_term=sum(P[i,j]*np.kron(varphis[i],np.conj(varphis[j]).T) for i in range(len(varphis)) for j in range(len(varphis)))
    lambda_=np.linalg.norm(first_term-second_term)
    return c_a,P,lambda_


def Optimization_final(Oas,varphis):
    number_oas = len(Oas)
    number_phis=len(varphis)
    dims=np.shape(Oas[0])[0]
    c_a=cp.Variable(number_oas)
    lambda_=cp.Variable()
    P=cp.Variable((number_phis,number_phis),complex=True)
    P.value=np.eye(number_phis)
    constraints=[P>>np.eye(number_phis)]
    first_term=sum(c_a[i]*Oas[i] for i in range(number_oas)) 
    second_term=sum(P[i,j]*np.kron(varphis[i],np.conj(varphis[j]).T) for i in range(number_phis) for j in range(number_phis))
    constraints+=[lambda_*np.eye(dims)>>first_term-second_term,first_term-second_term>>-lambda_*np.eye(dims)]
    objective=cp.Minimize(lambda_)
    prob=cp.Problem(objective,constraints)
    print("HAHA")
    prob.solve()
    print("The minimum of lambda is {}".format(prob.value))
    return c_a.value,P.value,lambda_.value


def construct_phis(M):
    # M is the mapping from boundary space to physical space. Here the column of M is the physical space and the row is the boundary space.
    orthogonal_complement=list(scipy.linalg.null_space(M.T).T)
    return_list=[]
    for vector in orthogonal_complement:
        return_list.append(vector.reshape(-1,1))
    return return_list


def D1MPS(group_num,MPS,D,chi):
    # The index of MPS is defined as 0:physics, 1:left,2:right
    start_MPS=MPS
    for i in range(1,group_num):
        start_MPS=ncon([start_MPS,MPS],[[-j for j in range(1,i+1)]+[-(i+2),1],[-(i+1),1,-(i+3)]])
    return start_MPS.reshape([D**group_num,chi**2])

def D2PEPS(H,V,D,chi,PEPS):
    # The physical index is from left to right, and bottom to top.
    horizontal_PEPS=PEPS
    for i in range(1,H):
        horizontal_PEPS=ncon([horizontal_PEPS,PEPS],[[-j for j in range(1,i+1)]+[-j for j in range(i+2,2*i+3)]+[1]+[-j for j in range(2*i+6,3*i+6)],[-(i+1),1,-(2*i+3),-(2*i+4),-(2*i+5)]])
    start_PEPS=deepcopy(horizontal_PEPS)
    for i in range(1,V):
        start_PEPS=ncon([start_PEPS,horizontal_PEPS],[[-j for j in range(1,H*i+1)]+[-j for j in range((i+1)*H+1,(i+1)*H+(i+1))]+[-j for j in range((i+1)*H+(i+2),(i+2)*H+(i+2))]+[-j for j in range((i+2)*H+(i+2),(i+2)*H+(2*i+2))]+[j for j in range(1,H+1)],[-j for j in range(H*i+1,H*(i+1)+1)]+[-((i+1)*H+(i+1))]+[H-j for j in range(H)]+[-((i+2)*H+(2*i+2))]+[-j for j in range((i+2)*H+(2*i+3),(i+3)*H+(2*i+3))]])
    return start_PEPS.reshape([D**(V*H),chi**(2*H+2*V)])


def testMPSExample(filename=None):
    MPS=np.zeros([3,2,2])
    sigmaminus=np.array([[0,0],[1,0]])
    sigmaplus=np.array([[0,1],[0,0]])
    sigmaz=np.array([[1,0],[0,-1]])
    MPS[0,:,:]=np.sqrt(2)*sigmaminus
    MPS[1,:,:]=-sigmaz
    MPS[2,:,:]=-np.sqrt(2)*sigmaplus
    M1=D1MPS(4,MPS,3,2)
    varphis=construct_phis(M1)
    Oas=[]
    import qutip as qt
    X=qt.spin_Jx(1)
    Y=qt.spin_Jy(1)
    Z=qt.spin_Jz(1)
    SS=np.kron(X,X)+np.kron(Y,Y)+np.kron(Z,Z)
    S1S2=np.kron(SS,np.eye(9))
    S2S3=np.kron(np.kron(np.eye(3),SS),np.eye(3))
    S3S4=np.kron(np.eye(9),SS)
    #S1S2=np.kron(SS,np.eye(3))
    #S2S3=np.kron(np.eye(3),SS)
    Oas.append(np.eye(81))
    #Oas.append(np.eye(27))
    Oas.append(S1S2)
    Oas.append(S2S3)
    Oas.append(S3S4)
    Oas.append(S1S2.dot(S1S2))
    Oas.append(S2S3.dot(S2S3))
    Oas.append(S3S4.dot(S3S4))
    c_a,P,lambda_=Optimization_by_hand_V2(Oas,varphis)
    print(c_a,flush=True)
    test_found_operator(c_a,Oas,varphis,M1)
    if filename is not None:
        np.savez(filename,c_a=c_a,P=P,lambda_=lambda_)


def testD2PEPS():
    MT=toriccode_MT()
    PEPS=statisticalPEPS_PEPS(MT,2)
    M1=D2PEPS(3,2,4,2,PEPS)
    print(np.shape(M1),flush=True)
    Oas=[]
    X=np.array([[0,1],[1,0]])
    Z=np.array([[1,0],[0,-1]])
    I=np.eye(2)
    Oas.append(np.eye(2**12))
    Oas.append(kron_product([Z,Z,Z,I,I,I,I,Z,I,I,I,I]))
    Oas.append(kron_product([I,I,Z,Z,Z,I,I,I,I,Z,I,I]))
    Oas.append(kron_product([I,I,X,I,I,I,I,X,X,X,I,I]))
    Oas.append(kron_product([I,I,I,I,X,I,I,I,I,X,X,X]))
    Oas.append(kron_product([X,X,Z,I,I,I,I,Z,I,I,I,I]))
    Oas.append(kron_product([I,I,X,X,Z,I,I,I,I,Z,I,I]))
    Oas.append(kron_product([I,I,Z,I,I,I,I,Z,X,X,I,I]))
    Oas.append(kron_product([I,I,I,I,Z,I,I,I,I,X,X,Z]))
    #Oas.append(kron_product([X,X,Z,I,I,Z,I,I]))
    #Oas.append(kron_product([I,I,X,I,I,X,X,X]))
    #Oas.append(kron_product([Z,Z,Z,I,I,Z,I,I]))
    #Oas.append(kron_product([I,I,Z,I,I,Z,X,X]))
    for i in Oas:
        print(np.linalg.norm((i-Oas[0]).dot(M1)),flush=True)


def delta_tensor(D):
    delta_tensor=np.zeros([D,D,D],dtype=complex)
    for i in range(D):
        delta_tensor[i,i,i]=1
    return delta_tensor

def statisticalPEPS_PEPS(MT,D):
    #MT is a matrix which defined in the overleaf
    MT=MT.reshape([D,D,D,D])
    DT=delta_tensor(D)
    return_tensor=ncon([DT,DT,MT],[[-1,-3,1],[-2,-4,2],[1,2,-5,-6]])
    return return_tensor.reshape([D**2,D,D,D,D])

def toriccode_MT():
    MT=np.zeros([2,2,2,2])
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    if (i+j+k+l)%2==0:
                        MT[i,j,k,l]=1
    return MT/np.sqrt(2)


def test_a_class_of_Dual_MT():
    x1=np.random.random()*np.sqrt(2)
    x2=np.random.random()*np.sqrt(2)
    x3=np.random.random()*np.sqrt(2)
    r=np.sqrt(2)/(x1*x2*x3)
    theta=2*np.pi*np.random.random()
    u4=r*np.cos(theta)
    lambda4=r*np.sin(theta)
    sign1=np.random.choice([-1,1])
    sign2=np.random.choice([-1,1])
    sign3=np.random.choice([-1,1])
    P1=np.array([[x1,0],[0,sign1*x1]])
    P2=np.array([[x2,0],[0,sign2*x2]])
    P3=np.array([[x3,0],[0,sign3*x3]])
    P4=np.array([[lambda4,0],[0,u4]])
    toricMT=toriccode_MT()
    new_MT=ncon([P1,P2,P3,P4,toricMT],[[-1,1],[-2,2],[-3,3],[-4,4],[1,2,3,4]])
    print(np.square(new_MT).reshape(4,4))
    Op_v=P2@P4
    Op_h=P1@P3
    P1I=np.linalg.inv(P1)
    P2I=np.linalg.inv(P2)
    P3I=np.linalg.inv(P3)
    P4I=np.linalg.inv(P4)
    Op_vI=np.linalg.inv(Op_v)
    Op_hI=np.linalg.inv(Op_h)
    final_matrix=specific_contraction_statistical_PEPS(new_MT)
    print(np.shape(final_matrix),flush=True)
    varphis=construct_phis(final_matrix)
    Oas=[]
    X=np.array([[0,1],[1,0]])
    Z=np.array([[1,0],[0,-1]])
    I=np.eye(2)   
    Oas.append(kron_product([P1@Z@P1I,P2@Z@P2I,Op_h@Z@Op_hI,I,I,I,Op_v@Z@Op_vI,I,I,I,I,I]))
    Oas.append(kron_product([I,I,Op_h@Z@Op_hI,P2@Z@P2I,P3@Z@P3I,I,I,I,I,Op_v@Z@Op_vI,I,I]))
    Oas.append(kron_product([I,I,I,I,I,P1@Z@P1I,Op_v@Z@Op_vI,Op_h@Z@Op_hI,P4@Z@P4I,I,I,I]))
    Oas.append(kron_product([I,I,I,I,I,I,I,Op_h@Z@Op_hI,I,Op_v@Z@Op_vI,P3@Z@P3I,P4@Z@P4I]))
    Oas.append(kron_product([I,I,Op_v@X@Op_vI,I,I,I,Op_h@X@Op_hI,Op_v@X@Op_vI,I,Op_h@X@Op_hI,I,I]))
    finding_op=2*np.eye(2**12)-0.25*Oas[0]-0.25*Oas[1]-0.25*Oas[2]-0.25*Oas[3]-Oas[4]
    print(np.linalg.norm(finding_op.dot(final_matrix)),flush=True)
    orthogonal_matrix=np.hstack(varphis)
    print(np.linalg.norm(finding_op.dot(orthogonal_matrix)),flush=True)
    print(scipy.linalg.eigvalsh(projection_to_orthogonal),flush=True)
    
    






def testPEPSExample(filename=None):
    MT=toriccode_MT()
    PEPS=statisticalPEPS_PEPS(MT,2)
    M1=D2PEPS(2,2,4,2,PEPS)
    print(np.shape(M1),flush=True)
    varphis=construct_phis(M1)
    Oas=[]
    X=np.array([[0,1],[1,0]])
    Z=np.array([[1,0],[0,-1]])
    Oas.append(np.eye(2**8))
    I=np.eye(2)
    #Oas.append(kron_product([Z,Z,Z,I,I,I,I,Z,I,I,I,I]))
    Oas.append(kron_product([Z,Z,Z,I,I,Z,I,I]))
    #Oas.append(kron_product([I,I,Z,Z,Z,I,I,I,I,Z,I,I]))
    #Oas.append(kron_product([I,I,X,I,I,I,I,X,X,X,I,I]))
    Oas.append(kron_product([I,I,X,I,I,X,X,X]))
   #Oas.append(kron_product([I,I,I,I,X,I,I,I,I,X,X,X]))
    #Oas.append(kron_product([X,X,Z,I,I,I,I,Z,I,I,I,I]))
    Oas.append(kron_product([X,X,Z,I,I,Z,I,I]))
    #Oas.append(kron_product([I,I,X,X,Z,I,I,I,I,Z,I,I]))
    #Oas.append(kron_product([I,I,Z,I,I,I,I,Z,X,X,I,I]))
    Oas.append(kron_product([I,I,Z,I,I,Z,X,X]))
    #Oas.append(kron_product([I,I,I,I,Z,I,I,I,I,X,X,Z]))
    c_a,P,lambda_=Optimization_final(Oas,varphis)
    if filename is not None:
        np.savez(filename,c_a=c_a,P=P,lambda_=lambda_)

def specific_contraction_statistical_PEPS(MT):
    #    8     11
    # 5     7     10
    #    6     9 
    # 0     2     4
    #    1     3
    MT=MT.reshape([2,2,2,2])
    DT=delta_tensor(2)
    PEPS1=ncon([DT,DT,MT],[[-1,-3,1],[-2,-4,2],[1,2,-5,-6]])#0,1
    PEPS2=ncon([DT,DT,DT,MT],[[-1,-4,1],[-2,-5,2],[-3,-6,3],[1,2,3,-7]])#2,3,4
    PEPS3=ncon([DT,DT,DT,DT,MT],[[-1,-5,1],[-2,-6,2],[-3,-7,3],[-4,-8,4],[1,2,3,4]])#5,6,7,8
    PEPS4=ncon([DT,DT,DT,MT],[[-1,-5,1],[-2,-6,2],[-3,-7,3],[-4,1,2,3]])#9,10,11
    final_PEPS=ncon([PEPS1,PEPS2,PEPS3,PEPS4],[[-1,-2,-13,-14,1,2],[-3,-4,-5,1,-15,-16,3],[-6,-7,-8,-9,-17,2,4,-18],[-10,-11,-12,4,3,-19,-20]])
    return final_PEPS.reshape([2**12,2**8])

def test_statistical_model(file_name=None):
    MT=toriccode_MT()
    final_matrix=specific_contraction_statistical_PEPS(MT)
    print(np.shape(final_matrix),flush=True)
    varphis=construct_phis(final_matrix)
    Oas=[]
    X=np.array([[0,1],[1,0]])
    Z=np.array([[1,0],[0,-1]])
    Oas.append(np.eye(2**12))
    I=np.eye(2)   
    Oas.append(kron_product([Z,Z,Z,I,I,I,Z,I,I,I,I,I]))
    Oas.append(kron_product([I,I,Z,Z,Z,I,I,I,I,Z,I,I]))
    Oas.append(kron_product([I,I,I,I,I,Z,Z,Z,Z,I,I,I]))
    Oas.append(kron_product([I,I,I,I,I,I,I,Z,I,Z,Z,Z]))
    Oas.append(kron_product([I,I,X,I,I,I,X,X,I,X,I,I]))
    Oas.append(kron_product([Z,X,X,I,I,I,Z,I,I,I,I,I]))
    Oas.append(kron_product([I,I,X,X,Z,I,I,I,I,Z,I,I]))
    Oas.append(kron_product([I,I,I,I,I,X,X,Z,Z,I,I,I]))
    Oas.append(kron_product([I,I,I,I,I,I,I,X,I,X,Z,Z]))
    Oas.append(kron_product([I,I,Z,I,I,I,Z,X,I,X,I,I]))
    c_a,P,lambda_=Optimization_by_hand_V2(Oas,varphis)
    print(c_a,flush=True)
    test_found_operator(c_a,Oas,varphis,final_matrix)
    if file_name is not None:
        np.savez(file_name,c_a=c_a,P=P,lambda_=lambda_)

def statistical_dual_ISOTNS_MT():
    [alpha,beta]=(np.random.random([2]))*2.0*np.pi
    thetas=np.random.random([8])*2*np.pi
    phis=np.random.random([16])*2*np.pi
    C_array=np.zeros([2]*8,dtype=complex)

    theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, theta_7, theta_8 = thetas
    phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, phi_8, phi_9, phi_10, phi_11, phi_12, phi_13, phi_14, phi_15, phi_16 = phis

    c1 = np.cos(alpha)*np.cos(theta_1)*np.exp(1j*phi_1)
    c2 = np.cos(alpha)*np.sin(theta_1)*np.exp(1j*phi_2)
    c3 = np.sin(alpha)*np.cos(theta_2)*np.exp(1j*phi_3)
    c4 = np.sin(alpha)*np.sin(theta_2)*np.exp(1j*phi_4)

    c5 = np.cos(beta)*np.cos(theta_5)*np.exp(1j*phi_9)
    c6 = np.cos(beta)*np.sin(theta_5)*np.exp(1j*phi_10)
    c7 = np.sin(beta)*np.cos(theta_6)*np.exp(1j*phi_11)
    c8 = np.sin(beta)*np.sin(theta_6)*np.exp(1j*phi_12)

    c9 = np.sin(alpha)*np.cos(theta_3)*np.exp(1j*phi_5)
    c10 = np.sin(alpha)*np.sin(theta_3)*np.exp(1j*phi_6)
    c11 = np.cos(alpha)*np.cos(theta_4)*np.exp(1j*phi_7)
    c12 = np.cos(alpha)*np.sin(theta_4)*np.exp(1j*phi_8)

    c13 = np.sin(beta)*np.cos(theta_7)*np.exp(1j*phi_13)
    c14 = np.sin(beta)*np.sin(theta_7)*np.exp(1j*phi_14)
    c15 = np.cos(beta)*np.cos(theta_8)*np.exp(1j*phi_15)
    c16 = np.cos(beta)*np.sin(theta_8)*np.exp(1j*phi_16)
    W = np.array([
        [c1, c2, c3, c4],
        [c5, c6, c7, c8],
        [c9, c10, c11, c12],
        [c13, c14, c15, c16]
    ])
    return W

def statistical_model_parent_Hamiltonian(filename=None):
    MT=statistical_dual_ISOTNS_MT()
    Oas=[]               
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    SU2matrixs=[X,Z,Y]
    I=np.eye(2)
    M1=specific_contraction_statistical_PEPS(MT)
    print("M1 completed",flush=True)
    print(np.shape(M1),flush=True)
    varphis=construct_phis(M1)
    print("varphis completed",flush=True)
    Oas.append(np.eye(2**12))
    #Oas.append(kron_product([Z,Z,Z,I,I,I,Z,I,I,I,I,I]))
    #Oas.append(kron_product([I,I,Z,Z,Z,I,I,I,I,Z,I,I]))
    #Oas.append(kron_product([I,I,I,I,I,Z,Z,Z,Z,I,I,I]))
    #Oas.append(kron_product([I,I,I,I,I,I,I,Z,I,Z,Z,Z]))
    #Oas.append(kron_product([I,I,X,I,I,I,X,X,I,X,I,I]))
    #Oas.append(kron_product([Z,X,X,I,I,I,Z,I,I,I,I,I]))
    #Oas.append(kron_product([I,I,X,X,Z,I,I,I,I,Z,I,I]))
    #Oas.append(kron_product([I,I,I,I,I,X,X,Z,Z,I,I,I]))
    #Oas.append(kron_product([I,I,I,I,I,I,I,X,I,X,Z,Z]))
    #Oas.append(kron_product([I,I,Z,I,I,I,Z,X,I,X,I,I]))
    #PEPS=statisticalPEPS_PEPS(MT,2)
    #M1=D2PEPS(3,2,4,2,PEPS)
    #varphis=construct_phis(M1)
    #Oas.append(np.eye(2**12))
    for i1 in range(2):
        for i2 in range(2):
            for i3 in range(2):
                for i4 in range(2):
                    Op1=SU2matrixs[i1]
                    Op2=SU2matrixs[i2]
                    Op3=SU2matrixs[i3]
                    Op4=SU2matrixs[i4]
                    temp_op=kron_product([Op1,Op2,Op3,I,I,I,Op4,I,I,I,I,I])
                    temp_op+=kron_product([I,I,Op1,Op2,Op3,I,I,I,I,Op4,I,I])
                    temp_op+=kron_product([I,I,I,I,I,Op1,Op2,Op3,Op4,I,I,I])
                    temp_op+=kron_product([I,I,I,I,I,I,I,Op1,I,Op2,Op3,Op4])
                    #Oas.append(temp_op)
                    #Oas.append(kron_product([I,I,Op1,I,I,I,Op2,Op3,I,Op4,I,I]))
    c_a,P,lambda_=Optimization_by_hand_V2(Oas,varphis)
    if filename is not None:
        np.savez(filename,c_a=c_a,P=P,lambda_=lambda_)
    test_found_operator(c_a,Oas,varphis,M1)

if __name__=="__main__":
    #test_statistical_model("./data/Toriccode_statistical_model.npz")
    #testMPSExample()
    statistical_model_parent_Hamiltonian("./data/dual_statistical_model_parent_Hamiltonian.npz")






    

    