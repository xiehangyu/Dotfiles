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
    project_channel=vectorized_channel
    project_channel -= kron_product([np.eye(chi**2),np.eye(chi**2),PO1,PO1]).dot(vectorized_channel)
    project_channel -= kron_product([PO1,np.eye(chi**2),np.eye(chi**2),PO1]).dot(vectorized_channel)
    project_channel += kron_product([PO1,np.eye(chi**2),PO1,PO1]).dot(vectorized_channel)
    project_channel += kron_product([state1]*4)
    return project_channel.reshape([chi**2]*4)

def check_condition_from_foldedchannel(fold_c,chi):
    right_c=fold_c.reshape([chi**4,chi**4])
    right_c=np.dot(right_c,right_c.conj().T)
    print("The right deviation is {}".format(LA.norm(right_c-np.eye(chi**4))))
    left_c=fold_c.transpose([0,3,2,1]).reshape([chi**4,chi**4])
    left_c=np.dot(left_c.conj().T,left_c)
    print("The left deviation is {}".format(LA.norm(left_c-np.eye(chi**4))))

    


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

