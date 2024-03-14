import numpy as np
from ncon import ncon 
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.stats import unitary_group
from scipy.linalg import expm
from scipy.optimize import fsolve
from scipy.optimize import approx_fprime
def R_derivative(C_array,l1,l2,d1,d2,alpha,l,d,r,u):
    if (l1 != l) or (d1 != d):
        return 0
    else:
        return np.conjugate(C_array[alpha,l2,d2,r,u])


def R_derivative_compl(C_array,l1,l2,d1,d2,alpha,l,d,r,u):
    if (l2 != l) or (d2 != d):
        return 0
    else:
        return C_array[alpha,l1,d1,r,u]

def R_derivative_re(C_array,l1,l2,d1,d2,alpha,l,d,r,u):
    return R_derivative(C_array,l1,l2,d1,d2,alpha,l,d,r,u) + R_derivative_compl(C_array,l1,l2,d1,d2,alpha,l,d,r,u)

def R_derivative_im(C_array,l1,l2,d1,d2,alpha,l,d,r,u):
    return 1j*(R_derivative(C_array,l1,l2,d1,d2,alpha,l,d,r,u) - R_derivative_compl(C_array,l1,l2,d1,d2,alpha,l,d,r,u))

def L_derivative(C_array,d1,d2,r1,r2,alpha,l,d,r,u):
    if (d1 != d) or (r1 != r):
        return 0
    else:
        return np.conjugate(C_array[alpha,l,d2,r2,u])

def L_derivative_compl(C_array,d1,d2,r1,r2,alpha,l,d,r,u):
    if (d2 != d) or (r2 != r):
        return 0
    else:
        return C_array[alpha,l,d1,r1,u]

def L_derivative_re(C_array,d1,d2,r1,r2,alpha,l,d,r,u):
    return L_derivative(C_array,d1,d2,r1,r2,alpha,l,d,r,u) + L_derivative_compl(C_array,d1,d2,r1,r2,alpha,l,d,r,u)

def L_derivative_im(C_array,d1,d2,r1,r2,alpha,l,d,r,u):
    return 1j*(L_derivative(C_array,d1,d2,r1,r2,alpha,l,d,r,u) - L_derivative_compl(C_array,d1,d2,r1,r2,alpha,l,d,r,u))

def CM_dualiso(C_array,D,chi):
    CM=np.zeros([2*chi**4,2*D*chi**4],dtype=complex)
    rowindex=0
    for l1 in range(chi):
        for d1 in range(chi):
            for l2 in range(chi):
                for d2 in range(chi):
                    if(l2*chi+d2<l1*chi+d1):
                        continue
                    colindex=0
                    for alpha in range(D):
                        for l in range(chi):
                            for d in range(chi):
                                for r in range(chi):
                                    for u in range(chi):
                                        CM[rowindex,colindex]=np.real(R_derivative_re(C_array,l1,l2,d1,d2,alpha,l,d,r,u))
                                        colindex+=1
                                        CM[rowindex,colindex]=np.real(R_derivative_im(C_array,l1,l2,d1,d2,alpha,l,d,r,u))
                                        colindex+=1
                    rowindex+=1
                    if(l2*chi+d2>l1*chi+d1):
                        colindex=0
                        for alpha in range(D):
                            for l in range(chi):
                                for d in range(chi):
                                    for r in range(chi):
                                        for u in range(chi):
                                            CM[rowindex,colindex]=np.imag(R_derivative_re(C_array,l1,l2,d1,d2,alpha,l,d,r,u))
                                            colindex+=1
                                            CM[rowindex,colindex]=np.imag(R_derivative_im(C_array,l1,l2,d1,d2,alpha,l,d,r,u))
                                            colindex+=1
                        rowindex+=1
    for d1 in range(chi):
        for r1 in range(chi):
            for d2 in range(chi):
                for r2 in range(chi):
                    if(d2*chi+r2<d1*chi+r1):
                        continue
                    colindex=0
                    for alpha in range(D):
                        for l in range(chi):
                            for d in range(chi):
                                for r in range(chi):
                                    for u in range(chi):
                                        CM[rowindex,colindex]=np.real(L_derivative_re(C_array,d1,d2,r1,r2,alpha,l,d,r,u))
                                        colindex+=1
                                        CM[rowindex,colindex]=np.real(L_derivative_im(C_array,d1,d2,r1,r2,alpha,l,d,r,u))
                                        colindex+=1
                    rowindex+=1
                    if(d2*chi+r2>d1*chi+r1):
                        colindex=0
                        for alpha in range(D):
                            for l in range(chi):
                                for d in range(chi):
                                    for r in range(chi):
                                        for u in range(chi):
                                            CM[rowindex,colindex]=np.imag(L_derivative_re(C_array,d1,d2,r1,r2,alpha,l,d,r,u))
                                            colindex+=1
                                            CM[rowindex,colindex]=np.imag(L_derivative_im(C_array,d1,d2,r1,r2,alpha,l,d,r,u))
                                            colindex+=1
                        rowindex+=1

    print("The shape of the condition matrix is ({},{})".format(rowindex,colindex))
    return CM


def tangent_space_dual_ISOTNS(c_array,D,chi):
    total_free_parameter=2*chi**4*D
    CM=CM_dualiso(c_array,D,chi)
    rank=np.linalg.matrix_rank(CM)
    print("The tangent space dimension is {}".format(total_free_parameter-rank))
    return total_free_parameter-rank


def construct_equal_dual_ISOTNS_from_gate(gate,D):
    C_array=np.zeros([D,D,D,D,D],dtype=complex)
    for alpha in range(D):
        for l in range(D):
            for d in range(D):
                for r in range(D):
                    for u in range(D):
                        C_array[alpha,l,d,r,u]=gate[alpha*D**2+r*D+u,l*D**2+d*D]
    return C_array

def check_condition(C_array,D,chi):
    leftmatrix=np.zeros([chi**2,chi**2],dtype=complex)
    for l1 in range(chi):
        for d1 in range(chi):
            for l2 in range(chi):
                for d2 in range(chi):
                    for alpha in range(D):
                        for r in range(chi):
                            for u in range(chi):
                                leftmatrix[l1*chi+d1,l2*chi+d2]+=np.conjugate(C_array[alpha,l2,d2,r,u])*C_array[alpha,l1,d1,r,u]
    print("The left variation from identity is {}".format(np.linalg.norm(leftmatrix-np.identity(chi**2))))
    rightmatrix=np.zeros([chi**2,chi**2],dtype=complex)
    for d1 in range(chi):
        for r1 in range(chi):
            for d2 in range(chi):
                for r2 in range(chi):
                    for alpha in range(D):
                        for l in range(chi):
                            for u in range(chi):
                                rightmatrix[d1*chi+r1,d2*chi+r2]+=np.conjugate(C_array[alpha,l,d2,r2,u])*C_array[alpha,l,d1,r1,u]
    print("The right variation from identity is {}".format(np.linalg.norm(rightmatrix-np.identity(chi**2))))

def dual_unitary_gates(D):
    SWAP=np.zeros([D**2,D**2],dtype=complex)
    for i in range(D):
        for j in range(D):
            SWAP[i*D+j,j*D+i]=1
    Z=np.diag([(D-1)/2-i for i in range(D)])
    ZZ=np.kron(Z,Z)
    u1=unitary_group.rvs(D)
    u2=unitary_group.rvs(D)
    u3=unitary_group.rvs(D)
    u4=unitary_group.rvs(D)
    j1=np.random.random()
    return np.kron(u1,u2)@expm(4j*j1*np.pi*ZZ)@SWAP@np.kron(u3,u4)
def dual_sequential_gates(D):
    u1=unitary_group.rvs(D**2)
    u2=dual_unitary_gates(D)
    return_u=np.kron(u2,np.eye(D))@np.kron(np.eye(D),u1)
    u1=unitary_group.rvs(D)
    u2=unitary_group.rvs(D)
    u3=unitary_group.rvs(D)
    u4=unitary_group.rvs(D)
    u5=unitary_group.rvs(D)
    u6=unitary_group.rvs(D)
    return np.kron(np.kron(u1,u2),u3)@return_u@np.kron(np.kron(u4,u5),u6)


def classical_partition_function_C_array():
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
    W=W.flatten().reshape([2,2,2,2])
    for l in range(2):
        for d in range(2):
            for r in range(2):
                for u in range(2):
                    C_array[l,d,r,u,l,d,r,u]=W[l,d,r,u]
    return C_array.reshape([16,2,2,2,2])
    


def controlled_dual_unitary_gates(D):
    returnmatrix=np.zeros([D**3,D**3],dtype=complex)
    for i in range(D):
        temp_matrix=dual_unitary_gates(D)
        returnmatrix += np.kron(temp_matrix,np.diag([1 if i==j else 0 for j in range(D)]))
    u1=unitary_group.rvs(D)
    u2=unitary_group.rvs(D)
    u3=unitary_group.rvs(D)
    u4=unitary_group.rvs(D)
    u5=unitary_group.rvs(D)
    u6=unitary_group.rvs(D)
    return np.kron(np.kron(u1,u2),u3)@returnmatrix@np.kron(np.kron(u4,u5),u6)



def Check_gothrough_condition(C_array,chi):
    state_left=ncon([np.conj(C_array),C_array],[[1,-1,-3,-5,2],[1,-2,-4,-6,2]])
    state_right=ncon([np.conj(C_array),C_array/(chi+0.0),np.eye(chi)],[[1,-1,2,-5,3],[1,-2,2,-6,3],[-3,-4]])
    return np.linalg.norm(state_left-state_right)

def Naive_generation_go_through(D,chi):
    U=unitary_group.rvs(D*chi).reshape([D,chi,D,chi])
    M=(np.random.random([D**2,D**2])+1j*np.random.random([D**2,D**2])).reshape([D,D,D,D])
    state=np.array([1]+[0]*(D-1))
    return ncon([U,M,state],[[-1,-5,-3,1],[1,-4,-2,2],[2]])

def test_Naive_generation_go_through(D):
    Check_gothrough_condition(Naive_generation_go_through(D),D)

def constraint_functions(chi):
    functions=[]
    for l1 in range(chi):
        for d1 in range(chi):
            for r1 in range(chi):
                for l2 in range(l1,chi):
                    for d2 in range(chi):
                        for r2 in range(chi):
                            if l1*chi**2+r1*chi+d1>l2*chi**2+r2*chi+d2:
                                continue
                            if d1!=d2:
                                def func(C_arrayflatten,D,chi,l1=l1,d1=d1,r1=r1,l2=l2,d2=d2,r2=r2):  
                                    C_arrays=(C_arrayflatten).reshape([2,D,chi,chi,chi,chi])
                                    C_array=C_arrays[0]+1j*C_arrays[1]
                                    temp_matrix1=(C_array[:,l1,d1,r1,:]).reshape([D,chi])
                                    temp_matrix2=((np.conj(C_array))[:,l2,d2,r2,:]).reshape([D,chi])
                                    return np.real(np.trace(temp_matrix1@(temp_matrix2.T)))
                                functions.append(func)
                                def func(C_arrayflatten,D,chi,l1=l1,d1=d1,r1=r1,l2=l2,d2=d2,r2=r2):
                                    C_arrays=(C_arrayflatten).reshape([2,D,chi,chi,chi,chi])
                                    C_array=C_arrays[0]+1j*C_arrays[1]
                                    temp_matrix1=(C_array[:,l1,d1,r1,:]).reshape([D,chi])
                                    temp_matrix2=((np.conj(C_array))[:,l2,d2,r2,:]).reshape([D,chi])
                                    return np.imag(np.trace(temp_matrix1@(temp_matrix2.T)))
                                functions.append(func)
                            if d1==d2:
                                def func(C_arrayflatten,D,chi,l1=l1,d1=d1,r1=r1,l2=l2,d2=d2,r2=r2):
                                    C_arrays=(C_arrayflatten).reshape([2,D,chi,chi,chi,chi])
                                    C_array=C_arrays[0]+1j*C_arrays[1]
                                    temp_matrix1=(C_array[:,l1,:,r1,:]).reshape([D,chi,chi])
                                    temp_matrix2=((np.conj(C_array))[:,l2,:,r2,:]).reshape([D,chi,chi])
                                    referencevalue=1/(chi+0.0)*ncon([temp_matrix1,temp_matrix2],[[1,2,3],[1,2,3]])
                                    temp_matrix1=(temp_matrix1[:,d1,:]).reshape([D,chi])
                                    temp_matrix2=(temp_matrix2[:,d2,:]).reshape([D,chi])
                                    return np.real(np.trace(temp_matrix1@(temp_matrix2.T))-referencevalue)
                                functions.append(func)
                                if r1!=r2 or l1!=l2:
                                    def func(C_arrayflatten,D,chi,l1=l1,d1=d1,r1=r1,l2=l2,d2=d2,r2=r2):
                                        C_arrays=(C_arrayflatten).reshape([2,D,chi,chi,chi,chi])
                                        C_array=C_arrays[0]+1j*C_arrays[1]
                                        temp_matrix1=(C_array[:,l1,:,r1,:]).reshape([D,chi,chi])
                                        temp_matrix2=((np.conj(C_array))[:,l2,:,r2,:]).reshape([D,chi,chi])
                                        referencevalue=1/(chi+0.0)*ncon([temp_matrix1,temp_matrix2],[[1,2,3],[1,2,3]])
                                        temp_matrix1=(temp_matrix1[:,d1,:]).reshape([D,chi])
                                        temp_matrix2=(temp_matrix2[:,d2,:]).reshape([D,chi])
                                        return np.imag(np.trace(temp_matrix1@(temp_matrix2.T))-referencevalue)
                                    functions.append(func)
    return functions



def calculate_tangent_space_dimension_go_through(C_array,D,chi):
    point=np.concatenate([np.real(C_array).flatten(),np.imag(C_array).flatten()])
    functions=constraint_functions(chi)
    sum=0
    for func in functions:
        sum+=np.abs(func(point,D,chi))
    print("The sum of the constraint functions is {}".format(sum))
    epsilon=(np.finfo(float).eps)**3
    num_functions=len(functions)
    num_variables=len(point)
    M=np.zeros([num_functions,num_variables])
    for i, func in enumerate(functions):
        def wrapped_func(x):
            return func(x,D,chi)
        M[i,:] = approx_fprime(point, wrapped_func, epsilon)
    print(np.shape(M))#//it should be chi^4(chi^2-1)*(2Dchi^4)
    s,v,d=np.linalg.svd(M)
    rank=np.sum(v>(1E-6*v[0]))
    print("The tangent space dimension is {}".format(2*chi**4*D-rank))
    return M


def go_through_tangent_space_naive_case(D):
    C_array=Naive_generation_go_through(D)
    return(calculate_tangent_space_dimension_go_through(C_array,D,D))




                            
                


def sequential_2qubitgates():
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    I=np.identity(2)
    XXI=np.kron(np.kron(X,X),I)
    YYI=np.kron(np.kron(Y,Y),I)
    ZZI=np.kron(np.kron(Z,Z),I)
    IXX=np.kron(np.kron(I,X),X)
    IYY=np.kron(np.kron(I,Y),Y)
    IZZ=np.kron(np.kron(I,Z),Z)
    Js=np.random.random(6)*2*np.pi
    Js[3]=0
    Js[2]=np.pi/4
    Js[4]=np.pi/4
    intm=expm(1j*Js[3]*IXX)@expm(1j*Js[4]*IYY)@expm(1j*Js[5]*IZZ)@expm(1j*Js[0]*XXI)@expm(1j*Js[1]*YYI)@expm(1j*Js[2]*ZZI)
    u1=unitary_group.rvs(2)
    u2=unitary_group.rvs(2)
    u3=unitary_group.rvs(2)
    u4=unitary_group.rvs(2)
    u5=unitary_group.rvs(2)
    return np.kron(np.kron(u1,u2),u3)@intm@np.kron(np.kron(u4,u5),I)

def test_sequential_2qubitgates():
    gate=sequential_2qubitgates()
    C_array=construct_equal_dual_ISOTNS_from_gate(gate,2)
    check_condition(C_array,2,2)
    tangent_space_dual_ISOTNS(C_array,2,2)

def test_controlled_dual_unitary_gates(D):
    gate=controlled_dual_unitary_gates(D)
    C_array=construct_equal_dual_ISOTNS_from_gate(gate,D)
    check_condition(C_array,D,D)
    tangent_space_dual_ISOTNS(C_array,D,D)
def test_dual_sequential_gates(D):
    gate=dual_sequential_gates(D)
    C_array=construct_equal_dual_ISOTNS_from_gate(gate,D)
    check_condition(C_array,D,D)
    tangent_space_dual_ISOTNS(C_array,D,D)
    

def translation_gates(D):
    intm=np.zeros([D**3,D**3],dtype=complex)
    for i in range(D):
        for j2 in range(D):
            for k in range(D):
                intm[i*D**2+j2*D+k,j2*D**2+k*D+i]=1
    u1=unitary_group.rvs(D)
    u2=unitary_group.rvs(D)
    u3=unitary_group.rvs(D)
    u4=unitary_group.rvs(D)
    u5=unitary_group.rvs(D)
    u6=unitary_group.rvs(D)
    return np.kron(np.kron(u1,u2),u3)@intm@np.kron(np.kron(u4,u5),u6)

def test_classical_partition_function_C_array():
    C_array=classical_partition_function_C_array()
    check_condition(C_array,D=16,chi=2)
    tangent_space_dual_ISOTNS(C_array,D=16,chi=2)

def test_translation_gates(D):
    gate=translation_gates(D)
    C_array=construct_equal_dual_ISOTNS_from_gate(gate,D)
    check_condition(C_array,D,D)
    tangent_space_dual_ISOTNS(C_array,D,D)

def equations_for_dual_TNS(vars,D,chi):
    C_array=vars.reshape([D,chi,chi,chi,chi])
    eqs=[]
    for l1 in range(chi):
        for d1 in range(chi):
            for l2 in range(chi):
                for d2 in range(chi):
                    if(l2*chi+d2<l1*chi+d1):
                        continue
                    tempsum=0
                    for alpha in range(D):
                        for r in range(chi):
                            for u in range(chi):
                                tempsum+=np.conjugate(C_array[alpha,l2,d2,r,u])*C_array[alpha,l1,d1,r,u]
                    if l1==l2 and d1==d2:
                        eqs.append(tempsum-1)
                    else:
                        eqs.append(tempsum)
    for d1 in range(chi):
        for r1 in range(chi):
            for d2 in range(chi):
                for r2 in range(chi):
                    if(d2*chi+r2<d1*chi+r1):
                        continue
                    tempsum=0
                    for alpha in range(D):
                        for l in range(chi):
                            for u in range(chi):
                                tempsum+=np.conjugate(C_array[alpha,l,d2,r2,u])*C_array[alpha,l,d1,r1,u]
                    if d1==d2 and r1==r2:
                        eqs.append(tempsum-1)
                    else:
                        eqs.append(tempsum)
    return np.array(eqs+[0j]*(chi**4*D-len(eqs)))


def random_guess_tangent_space(D,chi):
    init_guess=np.random.random([D,chi,chi,chi,chi])+1j*np.random.random([D,chi,chi,chi,chi])
    init_guess=init_guess.flatten()
    C_array=fsolve(equations_for_dual_TNS,init_guess,args=(D,chi)).reshape([D,chi,chi,chi,chi])
    check_condition(C_array,D,chi)
    tangent_space_dual_ISOTNS(C_array,D,chi)



def minimization_of_costfunction(vars,D,chi):
    vars_temp=vars.reshape([2,D,chi,chi,chi,chi])
    vars_tensor=vars_temp[0]+1j*vars_temp[1]
    tempmatrix1=vars_tensor.transpose([1,2,0,3,4]).reshape([chi*chi,D*chi*chi])
    costfunction1=tempmatrix1@np.conj(tempmatrix1).T-np.identity(chi*chi,dtype=complex)
    cost1=np.trace(costfunction1@np.conj(costfunction1).T)
    tempmatrix2=vars_tensor.transpose(2,3,0,1,4).reshape([chi*chi,D*chi*chi])
    costfunction2=tempmatrix2@np.conj(tempmatrix2).T-np.identity(chi*chi,dtype=complex)
    cost2=np.trace(costfunction2@np.conj(costfunction2).T)
    return np.sqrt(np.real(cost1+cost2))

def guess_by_minimization(D,chi):
    init_guess=np.random.random([2,D,chi,chi,chi,chi])
    init_guess=init_guess.flatten()
    result=minimize(minimization_of_costfunction,init_guess,args=(D,chi),method='BFGS')
    print(minimization_of_costfunction(result.x,D,chi))
    C_array_temp=result.x.reshape([2,D,chi,chi,chi,chi])
    C_array=C_array_temp[0]+1j*C_array_temp[1]

    check_condition(C_array,D,chi)
    tangent_space_dual_ISOTNS(C_array,D,chi)

def go_through_minimization_costfunction(vars,D,chi):
    vars_temp=vars.reshape([2,D,chi,chi,chi,chi])
    vars_tensor=vars_temp[0]+1j*vars_temp[1]
    return Check_gothrough_condition(vars_tensor,chi)**0.5
    

def Go_through_guess_by_minimization(D,chi):
    init_guess=np.random.random([2,D,chi,chi,chi,chi])
    init_guess=init_guess.flatten()
    result=minimize(go_through_minimization_costfunction,init_guess,args=(D,chi),method='BFGS')
    print(go_through_minimization_costfunction(result.x,D,chi))
    C_array_temp=result.x.reshape([2,D,chi,chi,chi,chi])
    C_array=C_array_temp[0]+1j*C_array_temp[1]
    return(calculate_tangent_space_dimension_go_through(C_array,D,chi))
    