import numpy as np
import matplotlib.pyplot as plt
from finding import translation_operator
from scipy.stats import unitary_group
from scipy.linalg import expm
from finding import tilde_U
from sympy import Matrix



def dual_unitary_gate(J):
    sigma_z=np.array([[1,0],[0,-1]])
    SWAP=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    temp_gate=np.kron(sigma_z,sigma_z)*J*1j
    return SWAP.dot(expm(temp_gate))

def self_random_matrix():
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    x,y,z=(np.random.rand(3)-0.5)/2
    H=x*X+y*Y+z*Z
    return expm(-1j*H)


def self_kron(op_list):
    temp_op=op_list[0]
    for i in range(1,len(op_list)):
        temp_op=np.kron(temp_op,op_list[i])
    return temp_op

def two_site_operator_with_random(U):
    u1=self_random_matrix()
    u2=self_random_matrix()
    return np.kron(u1,u2).dot(U)

def basic_floquet_operator_random(U1,U2,L):
    '''
    U is the unitary operator
    L is half the number of qubits
    '''
    temp_operator_list=[] 
    for i in range(L):
        temp_operator_list.append(two_site_operator_with_random(U1))
    temp_operator1=self_kron(temp_operator_list)
    temp_operator_list=[]
    for i in range(L):
        temp_operator_list.append(two_site_operator_with_random(U2))
    temp_operator2=self_kron(temp_operator_list)
    translation_operator1=translation_operator(2*L)
    temp_operator2=translation_operator1.dot(temp_operator2.dot(translation_operator1.transpose()))
    return temp_operator2.dot(temp_operator1)


def calculate_SFF(U1,U2,t,L):
    floquet_operator=basic_floquet_operator_random(U1,U2,L)
    temp_number=np.trace(np.linalg.matrix_power(floquet_operator,t))
    return np.abs(temp_number)**2



def ensemble_average(t,L,sample_size=1000):
    temp_sum=0
    CNOT=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    CZ=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])
    one_dual_unitary=dual_unitary_gate(0.9)
    for i in range(sample_size):
        temp_sum+=calculate_SFF(CNOT,CNOT,t,L)
    return temp_sum/sample_size

def operatorQK(operator_list,k):
    t=len(operator_list)
    temp_operator=np.zeros([4**t,4**t],dtype=complex)
    for i in range(1,1+t):
        phase=np.exp(1j*2*np.pi*k*i/t)
        temp_operator+=operator_list[i-1]*phase
    return temp_operator/t

def plot_ensemble_average(L):
    fs=open("./plot_data/ensemble_average_twoCNOT_L="+str(L)+".txt","w")
    for i in range(2,20):
        temp_number=ensemble_average(i,L)
        fs.write(str(i)+" "+str(temp_number)+"\n")
    fs.close()

def check_operator_norm(t):
    extra_factor_matrix=np.array([[1,0,0,1],[0,1,1,0],[0,1,1,0],[1,0,0,1]])
    extra_matrix=np.eye(1)
    for i in range(t):
        extra_matrix=np.kron(extra_matrix,extra_factor_matrix)
    QKs=[]
    TOs=[]
    T=translation_operator(2*t)
    T=T.dot(T)
    for i in range(1,t+1):
        temp_TO=np.eye(4**t)
        for j in range(i):
            temp_TO=temp_TO.dot(T)
        TOs.append(temp_TO)
    print(TOs[t-1])
    for k in range(t):
        QKs.append(operatorQK(TOs,k))
    for i in range(t):
        print(np.trace(QKs[i].dot(extra_matrix))/np.trace(QKs[i]))


def check_left_norm(t):
    extra_factor_matrix=np.array([[1,0,0,1],[0,1,1,0],[0,1,1,0],[1,0,0,1]])
    extra_matrix=np.eye(1)
    for i in range(t):
        extra_matrix=np.kron(extra_matrix,extra_factor_matrix)
    QKs=[]
    TOs=[]
    T=translation_operator(2*t)
    T=T.dot(T)
    for i in range(1,t+1):
        temp_TO=np.eye(4**t)
        for j in range(i):
            temp_TO=temp_TO.dot(T)
        TOs.append(temp_TO)
    print(TOs[t-1])
    for k in range(t):
        QKs.append(operatorQK(TOs,k))
    for i in range(t):
        normalization_factor=np.sqrt(np.trace(QKs[i]))
        print(normalization_factor)
        temp_matrix=QKs[i].dot(extra_matrix)
        temp_factor=np.trace(temp_matrix)/normalization_factor**2
        print(temp_factor)
        rest_matrix=temp_matrix/normalization_factor-temp_factor*QKs[i]/normalization_factor
        print(np.trace(rest_matrix.dot(rest_matrix.conjugate().transpose())))
    return QKs,extra_matrix

def generate_QKs(t):
    QKs=[]
    TOs=[]
    T=translation_operator(2*t)
    T=T.dot(T)
    for i in range(1,t+1):
        temp_TO=np.eye(4**t)
        for j in range(i):
            temp_TO=temp_TO.dot(T)
        TOs.append(temp_TO)
    for k in range(t):
        QKs.append(operatorQK(TOs,k))
    return QKs


def if_jordan(a):
    temp_matrix=Matrix(a)
    P,J=temp_matrix.jordan_form()
    j=np.array(J)
    temp_numpyarray=j-np.diag(np.diag(j))
    return np.sum(np.abs(temp_numpyarray))


def generate_project_operator_single(M, random_sample=20):
    '''
    M is the matrix
    this is acting on state
    '''
    temp_M=np.copy(M)
    projector=np.zeros(np.shape(temp_M))
    i=0
    for theta in np.linspace(0,2*np.pi,random_sample+1):
        if i==random_sample:
            break
        projector+=expm(temp_M*theta*1j)
        i=i+1
    projector=projector/random_sample
    return projector


def project_to_commutate_space(rho,M, random_sample=20):
    temp_M=np.copy(M)
    i=0
    final_state=np.zeros(np.shape(rho),dtype=complex)
    for theta in np.linspace(0,2*np.pi,random_sample+1):
        if i==random_sample:
            break
        projector=expm(temp_M*theta*1j)        
        final_state+=projector.dot(rho.dot(projector.conjugate().transpose()))
        i=i+1
    final_state=final_state/random_sample
    return final_state

def project_to_commutate_space_series(rho, Ms, repeated_num=5,random_sample=20):
    temp_state=np.copy(rho)
    for i in range(repeated_num):
        for M in Ms:
            temp_state=project_to_commutate_space(temp_state,M,random_sample)
    return temp_state 

def commutate(A,B):
    return A.dot(B)-B.dot(A)

def check_commutate(A,B):
    temp=commutate(A,B)
    return np.trace(temp.conj().transpose().dot(temp))/(np.trace(A.conj().transpose().dot(A))*np.trace(B.conj().transpose().dot(B)))

def out_of_space_overlap_matrixform(A,B):
    '''
    calculate the component of B that is out of A
    '''
    rest_matrix=B-np.trace(A.conj().transpose().dot(B))*A/np.trace(A.conj().transpose().dot(A))
    return np.trace(rest_matrix.conj().transpose().dot(rest_matrix))/np.trace(B.conj().transpose().dot(B))

def square_norm(A):
    return np.trace(A.conj().transpose().dot(A))

def check_coefficient_outspace(A,ops):
    '''
    ops is a list of operators
    '''
    temp_A=np.copy(A)
    for op in ops:
        temp_A=temp_A-np.trace(op.conj().transpose().dot(temp_A))*op/np.trace(op.conj().transpose().dot(op))
    return np.trace(temp_A.conj().transpose().dot(temp_A))/np.trace(A.conj().transpose().dot(A))
        

def general_check_overlap(t,O):
    CNOT=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    CNOT_tilde=tilde_U(CNOT)
    CNOT_tilde_seri=np.eye(1)
    for i in range(t):
        CNOT_tilde_seri=np.kron(CNOT_tilde_seri,CNOT_tilde)
    I=np.eye(2)
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    TO=translation_operator(2*t)
    CNOT_tilde_seri_TO=TO.dot(CNOT_tilde_seri.dot(TO.conj().transpose()))
    XI=np.kron(X,I)
    YI=np.kron(Y,I)
    ZI=np.kron(Z,I)
    odd_site_X=np.eye(1)
    for i in range(t):
        odd_site_X=np.kron(odd_site_X,XI)
    odd_site_Y=np.eye(1)
    for i in range(t):
        odd_site_Y=np.kron(odd_site_Y,YI)
    odd_site_Z=np.eye(1)
    for i in range(t):
        odd_site_Z=np.kron(odd_site_Z,ZI)
    even_site_X=TO.dot(odd_site_X.dot(TO.conj().transpose()))
    even_site_Y=TO.dot(odd_site_Y.dot(TO.conj().transpose()))
    even_site_Z=TO.dot(odd_site_Z.dot(TO.conj().transpose()))
    project_ops=[odd_site_X,odd_site_Y,odd_site_Z,even_site_X,even_site_Y,even_site_Z]
    #O=O/np.sqrt(np.trace(O.conj().transpose().dot(O)))
    print("the original norm is {}".format(square_norm(O)))
    after_action=CNOT_tilde_seri.dot(O.dot(CNOT_tilde_seri.conj().transpose()))
    print("the out of space component after first CNOT is{}\n".format(out_of_space_overlap_matrixform(O,after_action)))
    print("the norm after first CNOT is {}\n".format(square_norm(after_action)))
    after_project=project_to_commutate_space_series(after_action,project_ops,repeated_num=2)
    for M in project_ops:
        print(check_commutate(after_project,M))
    print(np.trace(after_action.conj().transpose().dot(after_project))/np.trace(after_action.conj().transpose().dot(after_action)))
    print("The out of space component after first projection is{}\n".format(out_of_space_overlap_matrixform(O,after_project)))
    print("the norm after first projection is {}\n".format(square_norm(after_project)))
    after_action2=CNOT_tilde_seri_TO.dot(after_project.dot(CNOT_tilde_seri_TO.conj().transpose()))
    print("the out of space component after second CNOT is{}\n".format(out_of_space_overlap_matrixform(O,after_action2)))
    print("the norm after second CNOT is {}\n".format(square_norm(after_action2)))
    after_project2=project_to_commutate_space_series(after_action2,project_ops,repeated_num=2)
    for M in project_ops:
        print(check_commutate(after_project2,M))
    print(np.trace(after_action2.conj().transpose().dot(after_project2))/np.trace(after_action2.conj().transpose().dot(after_action2)))
    print("The out of space component after second projection is{}\n".format(out_of_space_overlap_matrixform(O,after_project2)))
    print("the norm after second projection is {}\n".format(square_norm(after_project2)))
    print("\n\n\n")
    return after_project2    

def check_overlap_two_self(t,O1,O2,matrixO):
    O1_tilde=tilde_U(O1)
    O2_tilde=tilde_U(O2)
    O1_tilde_seri=np.eye(1)
    for i in range(t):
        O1_tilde_seri=np.kron(O1_tilde_seri,O1_tilde)
    O2_tilde_seri=np.eye(1)
    for i in range(t):
        O2_tilde_seri=np.kron(O2_tilde_seri,O2_tilde)
    I=np.eye(2)
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    TO=translation_operator(2*t)
    O2_tilde_seri=TO.dot(O2_tilde_seri.dot(TO.conj().transpose()))
    XI=np.kron(X,I)
    YI=np.kron(Y,I)
    ZI=np.kron(Z,I)
    odd_site_X=np.eye(1)
    for i in range(t):
        odd_site_X=np.kron(odd_site_X,XI)
    odd_site_Y=np.eye(1)
    for i in range(t):
        odd_site_Y=np.kron(odd_site_Y,YI)
    odd_site_Z=np.eye(1)
    for i in range(t):
        odd_site_Z=np.kron(odd_site_Z,ZI)
    even_site_X=TO.dot(odd_site_X.dot(TO.conj().transpose()))
    even_site_Y=TO.dot(odd_site_Y.dot(TO.conj().transpose()))
    even_site_Z=TO.dot(odd_site_Z.dot(TO.conj().transpose()))
    project_ops=[odd_site_X,odd_site_Y,odd_site_Z,even_site_X,even_site_Y,even_site_Z]
    #O=O/np.sqrt(np.trace(O.conj().transpose().dot(O)))
    print("the original norm is {}".format(square_norm(matrixO)))
    after_action=O1_tilde_seri.dot(matrixO.dot(O1_tilde_seri.conj().transpose()))
    print("the out of space component after first gate is{}\n".format(out_of_space_overlap_matrixform(matrixO,after_action)))
    print("the norm after first gate is {}\n".format(square_norm(after_action)))
    after_project=project_to_commutate_space_series(after_action,project_ops,repeated_num=2)
    for M in project_ops:
        print(check_commutate(after_project,M))
    print(np.trace(after_action.conj().transpose().dot(after_project))/np.trace(after_action.conj().transpose().dot(after_action)))
    print("The out of space component after first projection is{}\n".format(out_of_space_overlap_matrixform(matrixO,after_project)))
    print("the norm after first projection is {}\n".format(square_norm(after_project)))
    after_action2=O2_tilde_seri.dot(after_project.dot(O2_tilde_seri.conj().transpose()))
    print("the out of space component after second gate is{}\n".format(out_of_space_overlap_matrixform(matrixO,after_action2)))
    print("the norm after second gate is {}\n".format(square_norm(after_action2)))
    after_project2=project_to_commutate_space_series(after_action2,project_ops,repeated_num=2)
    for M in project_ops:
        print(check_commutate(after_project2,M))
    print(np.trace(after_action2.conj().transpose().dot(after_project2))/np.trace(after_action2.conj().transpose().dot(after_action2)))
    print("The out of space component after second projection is{}\n".format(out_of_space_overlap_matrixform(matrixO,after_project2)))
    print("the norm after second projection is {}\n".format(square_norm(after_project2)))
    print("\n\n\n")
    return after_project2    


def check_overlap(t):
    QKs=generate_QKs(t)
    CNOT=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    CNOT_tilde=tilde_U(CNOT)
    CNOT_tilde_seri=np.eye(1)
    for i in range(t):
        CNOT_tilde_seri=np.kron(CNOT_tilde_seri,CNOT_tilde)
    I=np.eye(2)
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    TO=translation_operator(2*t)
    CNOT_tilde_seri_TO=TO.dot(CNOT_tilde_seri.dot(TO.conj().transpose()))
    extra_factor_matrix=np.array([[1,0,0,1],[0,1,1,0],[0,1,1,0],[1,0,0,1]])
    extra_matrix=np.eye(1)
    for i in range(t):
        extra_matrix=np.kron(extra_matrix,extra_factor_matrix)
    TO_extra_matrix=TO.dot(extra_matrix.dot(TO.conj().transpose()))
    XI=np.kron(X,I)
    YI=np.kron(Y,I)
    ZI=np.kron(Z,I)
    odd_site_X=np.eye(1)
    for i in range(t):
        odd_site_X=np.kron(odd_site_X,XI)
    odd_site_Y=np.eye(1)
    for i in range(t):
        odd_site_Y=np.kron(odd_site_Y,YI)
    odd_site_Z=np.eye(1)
    for i in range(t):
        odd_site_Z=np.kron(odd_site_Z,ZI)
    even_site_X=TO.dot(odd_site_X.dot(TO.conj().transpose()))
    even_site_Y=TO.dot(odd_site_Y.dot(TO.conj().transpose()))
    even_site_Z=TO.dot(odd_site_Z.dot(TO.conj().transpose()))
    project_ops=[odd_site_X,odd_site_Y,odd_site_Z,even_site_X,even_site_Y,even_site_Z]
    for i in range(t):
        '''
        print("now begin {}".format(i))
        QKs[i]=QKs[i]/np.sqrt(np.trace(QKs[i]))
        print("the original norm is {}".format(square_norm(QKs[i])))
        after_action=CNOT_tilde_seri.dot(QKs[i].dot(CNOT_tilde_seri.conj().transpose()))
        print("the out of space component after first CNOT is{}\n".format(out_of_space_overlap_matrixform(QKs[i],after_action)))
        print("the norm after first CNOT is {}\n".format(square_norm(after_action)))
        after_project=project_to_commutate_space_series(after_action,project_ops,repeated_num=2)
        for M in project_ops:
            print(check_commutate(after_project,M))
        print(np.trace(after_action.conj().transpose().dot(after_project))/np.trace(after_action.conj().transpose().dot(after_action)))
        print("The out of space component after first projection is{}\n".format(out_of_space_overlap_matrixform(QKs[i],after_project)))
        print("the norm after first projection is {}\n".format(square_norm(after_project)))
        after_action2=CNOT_tilde_seri_TO.dot(after_project.dot(CNOT_tilde_seri_TO.conj().transpose()))
        print("the out of space component after second CNOT is{}\n".format(out_of_space_overlap_matrixform(QKs[i],after_action2)))
        print("the norm after second CNOT is {}\n".format(square_norm(after_action2)))
        after_project2=project_to_commutate_space_series(after_action2,project_ops,repeated_num=2)
        for M in project_ops:
            print(check_commutate(after_project2,M))
        print(np.trace(after_action2.conj().transpose().dot(after_project2))/np.trace(after_action2.conj().transpose().dot(after_action2)))
        print("The out of space component after second projection is{}\n".format(out_of_space_overlap_matrixform(QKs[i],after_project2)))
        print("the norm after second projection is {}\n".format(square_norm(after_project2)))
        print("\n\n\n")
        '''
        print("now begin {}".format(i))
        QKs[i]=QKs[i]/np.sqrt(np.trace(QKs[i]))
        temp_t1=general_check_overlap(t,QKs[i])
        print("the second layer of {}".format(i))
        temp_t2=general_check_overlap(t,temp_t1)


CNOT=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])