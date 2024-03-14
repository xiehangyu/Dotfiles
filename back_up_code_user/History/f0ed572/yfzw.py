from calculate_dual_ISOTNS import *
import numpy.linalg as LA
def MaxEnS(chi):
    vectors = np.zeros([chi**2,1], dtype=complex)
    for i in range(chi):
        vectors[i*chi+i,0] = 1.0
    return vectors/np.sqrt(chi)

def vectorizedchannel(C_array,D,chi):
    C_matrix=np.reshape(C_array,[D,chi**4])
    channel=C_matrix.conj().T.dot(C_matrix)
    vectorized_channel=np.reshape(channel,[chi]*8)
    return vectorized_channel.transpose([0,4,1,5,2,6,3,7]).flatten()

def expect(OP,state):
    temp=np.real(np.vdot(state,OP.dot(state)))
    if temp < -1e-10:
        print(temp)
    if temp<=0 and temp > -1e-10:
        return 0
    else:
        return temp

def kron_product(a):
    result=a[0]
    for array in a[1:]:
        result=np.kron(result,array)
    return result

def project_to_dual(vectorized_channel,chi):
    PO1=MaxEnS(chi)
    state1=PO1
    PO1=PO1.dot(PO1.conj().T)
    project_channel=np.copy(vectorized_channel)
    project_channel.dtype=complex
    project_channel -= kron_product([np.eye(chi**2),np.eye(chi**2),PO1,PO1]).dot(vectorized_channel)
    project_channel -= kron_product([PO1,np.eye(chi**2),np.eye(chi**2),PO1]).dot(vectorized_channel)
    project_channel += kron_product([PO1,np.eye(chi**2),PO1,PO1]).dot(vectorized_channel)
    project_channel += kron_product([PO1,PO1,PO1,PO1]).dot(vectorized_channel)
    coef=np.vdot(kron_product([state1]*4),project_channel)
    if abs(coef)<1e-8:
        return False
    else:
        return project_channel.reshape([chi**2]*4)/coef


def random_PEPS(D,chi):
    return np.random.random([D,chi,chi,chi,chi])+1j*np.random.random([D,chi,chi,chi,chi])  

def random_dual_PEPS_fold(D,chi):
    vectorized_channel=vectorizedchannel(random_PEPS(D,chi),D,chi)
    fold_c=project_to_dual(vectorized_channel,chi)
    if fold_c is False:
        return False
    else:
        return fold_c

def random_dual_PEPS(D,chi):
    while True:
        PEPS=random_PEPS(D,chi)
        vectorized_channel=vectorizedchannel(PEPS,D,chi)
        fold_c=project_to_dual(vectorized_channel,chi)
        matrix_c=fold_c.reshape([chi]*8).transpose([0,2,4,6,1,3,5,7]).reshape([chi**4]*2)
        matrix_c=(matrix_c+matrix_c.conj().T)/2.0
        eigenvalues,eigenvectors=LA.eigh(matrix_c)
        #set negative values in eigenvalues to zero
        eigenvalues=np.real(eigenvalues)
        eigenvalues[eigenvalues<0]=0
        if (eigenvalues[-D]/eigenvalues[-1])<=1e-7:
            C_array=np.diag(np.sqrt(eigenvalues[-D:]))@(eigenvectors[:,-D:].conj().T)
            check_condition(C_array,D,chi)
            return C_array.reshape([D,chi,chi,chi,chi])
        else:
            print(eigenvalues[-D]/eigenvalues[-1])
            


def check_random_dual_PEPS_fold(D,chi):
    from contraction import check_condition_from_foldedchannel
    fold_c=random_dual_PEPS_fold(D,chi)
    if fold_c is False:
        print("The random PEPS channel cannot be dualized")
    else:
        check_condition_from_foldedchannel(fold_c,chi)
        print(decomposition_coefficient(fold_c.flatten(),chi))


def decomposition_coefficient(vectorized_channel,chi):
    PO1=MaxEnS(chi)
    PO1=PO1.dot(PO1.conj().T)
    PO2=np.eye(chi**2)-PO1
    POs=[PO1,PO2]
    decomposition_coefficient=np.zeros([2]*4)
    for l in range(2):
        POl=POs[l]
        for d in range(2):
            POd=POs[d]
            for r in range(2):
                POr=POs[r]
                for u in range(2):
                    POu=POs[u]
                    PO=np.kron(POl,np.kron(POd,np.kron(POr,POu)))
                    decomposition_coefficient[l,d,r,u]=np.sqrt(expect(PO,vectorized_channel))
    return decomposition_coefficient

def decomposition_coefficient_C_array(C_array,D,chi):
    return decomposition_coefficient(vectorizedchannel(C_array,D,chi),chi)


def distance1(C_array,D,chi):
    coef=decomposition_coefficient_C_array(C_array,D,chi)
    basis_value=np.abs(coef[0,0,0,0])
    return np.sqrt(LA.norm(coef)**2-basis_value**2)/basis_value


def check_decomposition(D,chi):
    init_guess=np.random.random([2,D,chi,chi,chi,chi])
    init_guess=init_guess.flatten()
    result=minimize(minimization_of_costfunction,init_guess,args=(D,chi),method='BFGS')
    C_array_temp=result.x.reshape([2,D,chi,chi,chi,chi])
    C_array=C_array_temp[0]+1j*C_array_temp[1]
    check_condition(C_array,D,chi)
    print(decomposition_coefficient_C_array(C_array,D,chi))


def distance1_sequential():
    gate=sequential_2qubitgates()
    C_array=construct_equal_dual_ISOTNS_from_gate(gate,2)
    check_condition(C_array,2,2)
    print(decomposition_coefficient_C_array(C_array,2,2))
    print("The distance is {}".format(distance1(C_array,2,2)))

def distance1_control_dual(D):
    gate=controlled_dual_unitary_gates(D)
    C_array=construct_equal_dual_ISOTNS_from_gate(gate,D)
    check_condition(C_array,D,D)
    print(decomposition_coefficient_C_array(C_array,D,D))
    print("The distance is {}".format(distance1(C_array,D,D)))

def distance1_permutation(D):
    gate=translation_gates(D)
    C_array=construct_equal_dual_ISOTNS_from_gate(gate,D)
    check_condition(C_array,D,D)
    print(decomposition_coefficient_C_array(C_array,D,D))
    print("The distance is {}".format(distance1(C_array,D,D)))

