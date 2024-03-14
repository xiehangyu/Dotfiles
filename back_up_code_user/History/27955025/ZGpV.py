import numpy as np
import scipy
from ncon import ncon
from calculate_dual_ISOTNS import *
import numpy.linalg as LA
from Moses_Move import *
from copy import deepcopy
def foldedchannel(C_array,D,chi):
    #return folded PEPS channel (fold_c)
    temp_matrix=C_array.reshape([D,chi**4])
    temp_matrix=temp_matrix.conj().T@temp_matrix
    return temp_matrix.reshape([chi]*8).transpose([0,4,1,5,2,6,3,7]).reshape([chi**2]*4)

def folded_operator_tensor(C_array,Op,D,chi):
    #return folded operator tensor (fold_o)
    temp_matrix=C_array.reshape([D,chi**4])
    temp_matrix=temp_matrix.conj().T@Op@temp_matrix
    return temp_matrix.reshape([chi]*8).transpose([0,4,1,5,2,6,3,7]).reshape([chi**2]*4)

def folded_identity(chi):
    #return folded identity (folded_in)
    return np.eye(chi).reshape([chi**2])/np.sqrt(chi)

def check_condition_from_foldedchannel(fold_c,chi):
    fold_in=folded_identity(chi)
    right=ncon([fold_c,fold_in,fold_in],[[-1,-2,1,2],[1],[2]])
    print("The right deviation is {}".format(LA.norm(right-ncon([fold_in,fold_in],[[-1],[-2]]))))
    left=ncon([fold_c,fold_in,fold_in],[[1,-1,-2,2],[1],[2]])
    print("The left deviation is {}".format(LA.norm(left-ncon([fold_in,fold_in],[[-1],[-2]]))))


def vertical_transform_matrix(fold_c,v_n,chi):
    # return vertical_transformation_matrix (ver_tran)
    folded_in=folded_identity(chi)
    new_c=ncon([fold_c,folded_in],[[-1,1,-2,-3],[1]])
    for i in range(v_n-1):
        new_c=ncon([fold_c,new_c],[[-1,1,-4-2*i,-5-2*i],[-2-j for j in range(i+1)]+[-3-i-j for j in range(i+1)]+[1]])
    new_c=ncon([new_c,folded_in],[[-j for j in range(1,2*v_n+1)]+[1],[1]])
    return new_c

def contraction_with_fixed_point(fold_c,Lx,chi):
    pass    

def horizontal_transform_matrix(fold_c,Lx,chi):
    #return horizontal_transformation_matrix (hor_tran)
    fold_in=folded_identity(chi)
    new_c=ncon([fold_c,fold_in],[[1,-1,-2,-3],[1]])
    for i in range(1,Lx):
        new_c=ncon([new_c,fold_c],[[-j for j in range(1,i+1)]+[1]+[-i-4-j for j in range(0,i)],[1,-i-1,-i-2,-i-3]])
    new_c=ncon([new_c,fold_in],[[-j for j in range(1,Lx+1)]+[1]+[-Lx-j for j in range(1,Lx+1)],[1]])
    return new_c

def left_boundary(ver_tran,v_n,chi):
    #return the transformation matrix contracted with the left identity boundary (left_c)
    folded_in=folded_identity(chi)
    return ncon([ver_tran]+[folded_in]*v_n,[[i for i in range(1,v_n+1)]+[-i for i in range(1,v_n+1)]]+[[i] for i in range(1,v_n+1)])


def contration_to_right(left_c,ver_tran,v_n):
    return ncon([left_c,ver_tran],[[i for i in range(1,v_n+1)],[v_n+1-i for i in range(1,v_n+1)]+[-i for i in range(1,v_n+1)]])

def right_boundary(left_c,v_n,chi):
    #Return the value (a number) of the result
    folded_in=folded_identity(chi)
    return ncon([left_c]+[folded_in]*v_n,[[i for i in range(1,v_n+1)]]+[[i] for i in range(1,v_n+1)])
        
def vertical_transform_with_operator(fold_c,fold_o,v_n,chi,posi):
    #here the position is calculated from 0 to v_n-1
    folded_in=folded_identity(chi)
    if posi==0:
        new_c=ncon([fold_o,folded_in],[[-1,1,-2,-3],[1]])
    else:
        new_c=ncon([fold_c,folded_in],[[-1,1,-2,-3],[1]])
    for i in range(v_n-1):
        if i==posi-1:
            new_c=ncon([fold_o,new_c],[[-1,1,-4-2*i,-5-2*i],[-2-j for j in range(i+1)]+[-3-i-j for j in range(i+1)]+[1]])
        else:
            new_c=ncon([fold_c,new_c],[[-1,1,-4-2*i,-5-2*i],[-2-j for j in range(i+1)]+[-3-i-j for j in range(i+1)]+[1]])
    new_c=ncon([new_c,folded_in],[[-j for j in range(1,2*v_n+1)]+[1],[1]])
    return new_c

def vertical_transform_with_operator_samecolumn(fold_c,fold_o1,fold_o2,v_n,chi,posi1,posi2):
    fold_in=folded_identity(chi)
    if posi1==0:
        new_c=ncon([fold_o1,fold_in],[[-1,1,-2,-3],[1]])
    elif posi2==0:
        new_c=ncon([fold_o2,fold_in],[[-1,1,-2,-3],[1]])
    else:
        new_c=ncon([fold_c,fold_in],[[-1,1,-2,-3],[1]])
    for i in range(v_n-1):
        if i==posi1-1:
            new_c=ncon([fold_o1,new_c],[[-1,1,-4-2*i,-5-2*i],[-2-j for j in range(i+1)]+[-3-i-j for j in range(i+1)]+[1]])
        elif i==posi2-1:
            new_c=ncon([fold_o2,new_c],[[-1,1,-4-2*i,-5-2*i],[-2-j for j in range(i+1)]+[-3-i-j for j in range(i+1)]+[1]])
        else:
            new_c=ncon([fold_c,new_c],[[-1,1,-4-2*i,-5-2*i],[-2-j for j in range(i+1)]+[-3-i-j for j in range(i+1)]+[1]])
    new_c=ncon([new_c,fold_in],[[-j for j in range(1,2*v_n+1)]+[1],[1]])
    return new_c
    

def PEPS_direct_contraction(fold_o1,fold_o2,x1,y1,x2,y2,fold_c,chi,Lx,Ly):
    #x1,x2,y1,y2 are all calculated from 0
    #Assume x2!=x1
    if x1!=x2:
        ver_tran=vertical_transform_matrix(fold_c,Ly,chi)
        ver_o1=vertical_transform_with_operator(fold_c,fold_o1,Ly,chi,y1)
        ver_o2=vertical_transform_with_operator(fold_c,fold_o2,Ly,chi,y2)
        if x1==0:
            start_c=ver_o1
        elif x2==0:
            start_c=ver_o2
        else:
            start_c=ver_tran
        left_c=left_boundary(start_c,Ly,chi)
        for i in range(1,Lx):
            if i==x1:
                left_c=contration_to_right(left_c,ver_o1,Ly)
            elif i==x2:
                left_c=contration_to_right(left_c,ver_o2,Ly)
            else:
                left_c=contration_to_right(left_c,ver_tran,Ly)
        return right_boundary(left_c,Ly,chi)
    else:
        ver_o1o2=vertical_transform_with_operator_samecolumn(fold_c,fold_o1,fold_o2,Ly,chi,y1,y2)
        ver_tran=vertical_transform_matrix(fold_c,Ly,chi)
        if x1==0:
            start_c=ver_o1o2
        else:
            start_c=ver_tran
        left_c=left_boundary(start_c,Ly,chi)
        for i in range(1,Lx):
            if i==x1:
                left_c=contration_to_right(left_c,ver_o1o2,Ly)
            else:
                left_c=contration_to_right(left_c,ver_tran,Ly)
        return right_boundary(left_c,Ly,chi)

def PEPS_direct_from_op(O1,O2,x1,y1,x2,y2,C_array,D,chi,Lx,Ly):
    #Assume x2>x1, x1,x2,y1,y2 are all calculated from 0
    fold_o1=folded_operator_tensor(C_array,O1,D,chi)
    fold_o2=folded_operator_tensor(C_array,O2,D,chi)
    fold_c=foldedchannel(C_array,D,chi)
    return PEPS_direct_contraction(fold_o1,fold_o2,x1,y1,x2,y2,fold_c,chi,Lx,Ly)


def MPS_fixedpoint_contraction(fold_o1,fold_o2,x2,y2,fold_c,chi):
    #Assume x2>0,y2>0,x1=y1=0
    Lx=x2+1
    htm=horizontal_transform_matrix(fold_c,Lx,chi)
    bd=chi**2
    htm=htm.transpose([j for j in range(Lx)]+[2*Lx-j-1 for j in range(Lx)])
    htm=htm.reshape([bd**Lx]*2).transpose()
    eigenvalues,eigenstates=scipy.linalg.eig(htm)
    idx=eigenvalues.argsort()[::-1]
    eigenvalues=eigenvalues[idx]
    eigenstates=eigenstates[:,idx]
    print(eigenvalues[:5])
    eigenstate=eigenstates[:,0]
    eigenstate=eigenstate.reshape([bd]*Lx)
    fold_in=folded_identity(chi)
    start_c=ncon([fold_o1,fold_in,fold_in],[1,-1,-2,2],[1,2])
    if x2>0:
        for i in range(1,x2):
            start_c=ncon([start_c,fold_c,fold_in],[[-j for j in range(1,i+1)]+[1],[1,-i-1,-i-2,2],[2]])
        if y2==0:
            start_c=ncon([start_c,fold_o2,fold_in,fold_in],[[-j for j in range(1,x2+1)]+[1],[1,-x2-1,2,3],[2],[3]])
            return ncon([start_c,eigenstate],[[i for i in range(1,Lx+1)],[j for j in range(Lx+1)]])
        start_c=ncon([start_c,fold_c,fold_in],[[-j for j in range(1,x2+1)]+[1],[1,-x2-1,2,-x2-2],[2]])
        start_c=ncon([start_c,eigenstate],[[i for i in range(1,Lx+1)]+[-1],[j for j in range(Lx+1)]])
        for i in range(y2-1):
            start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[2,1,3,-1],[2],[3]])
        return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[2,1,3,4],[2],[3],[4]])
    if x2==0:
        start_c=ncon([fold_o1,fold_in,eigenstate,fold_in],[[1,2,3,-1],[1],[2],[3]])
        for i in range(y2-1):
            start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[2,1,3,-1],[2],[3]])
        return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[2,1,3,4],[2],[3],[4]])
        


def MPS_direct_contraction(fold_o1,fold_o2,x2,y2,fold_c,chi):
    #Assume x1=y1=0,y2>=0, x2,y2 are both calculated from 0
    fold_in=folded_identity(chi)
    if x2>0:
        #Assume x1=y1=0, also assume x2>0, y2>=0, x2,y2 are both calculated from 0
        start_c=ncon([fold_o1,fold_in,fold_in,fold_in],[[1,2,-1,3],[1],[2],[3]])
        for i in range(x2-1):
            start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[1,2,-1,3],[2],[3]])
        if y2==0:
            return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[1,2,3,4],[2],[3],[4]])
        start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[1,2,3,-1],[2],[3]])
        for i in range(y2-1):
            start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[2,1,3,-1],[2],[3]])
        return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[2,1,3,4],[2],[3],[4]])
    if x2==0:
        start_c=ncon([fold_o1,fold_in,fold_in,fold_in],[[1,2,3,-1],[1],[2],[3]])
        for i in range(y2-1):
            start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[2,1,3,-1],[2],[3]])
        return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[2,1,3,4],[2],[3],[4]])
    if x2<0:
        x2=-x2
        start_c=ncon([fold_o1,fold_in,fold_in,fold_in],[[-1,1,2,3],[1],[2],[3]])
    for i in range(x2-1):
        start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[-1,2,1,3],[2],[3]])
    if y2==0:
        return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[2,3,1,4],[2],[3],[4]])
    start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[2,3,1,-1],[2],[3]])
    for i in range(y2-1):
        start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[2,1,3,-1],[2],[3]])
    return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[2,1,3,4],[2],[3],[4]])


def MPS_direct_from_op(O1,O2,x2,y2,C_array,D,chi):
    #Assume x1=y1=0, also assume x2>0, y2>=0, x2,y2 are both calculated from 0
    fold_o1=folded_operator_tensor(C_array,O1,D,chi)
    fold_o2=folded_operator_tensor(C_array,O2,D,chi)
    fold_c=foldedchannel(C_array,D,chi)
    return MPS_direct_contraction(fold_o1,fold_o2,x2,y2,fold_c,chi)

def general_MPS_direct_contraction_from_op(O1,O2,x1,y1,x2,y2,C_array,D,chi):
    #We assume y2>=y1
    if y1==0:
        return MPS_direct_from_op(O1,O2,x2-x1,y2-y1,C_array,D,chi)
    #We further assume y2>=y1 and x2>x1
    '''The code below if for Moses Move
    elif y2>=y1 and x2>x1:
        Lx=x2-x1+1
        iteration_y=y1
        based_C_arrays=[C_array]*Lx
        layer_to_moves=[C_array]*Lx
        for i in range(iteration_y):
            bottom_layer_C_arrays,upper_layer_C_arrays=Moses_Move(layer_to_moves,chi=chi,D=D)
            layer_to_moves=truncation_of_two_layers(upper_layer_C_arrays,deepcopy(based_C_arrays),chi=chi,D=D)
        fold_o1=folded_operator_tensor(layer_to_moves[0],O1,D,chi)
        fold_in=folded_identity(chi)
        start_c=ncon([fold_o1,fold_in,fold_in,fold_in],[[1,2,-1,3],[1],[2],[3]])
        for i in range(1,Lx-1):
            fold_c=foldedchannel(layer_to_moves[i],D,chi)
            start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[1,2,-1,3],[2],[3]])
        y=y2-y1
        if y==0:
            fold_o2=folded_operator_tensor(layer_to_moves[Lx-1],O2,D,chi)
            return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[1,2,3,4],[2],[3],[4]])
        else:
            fold_o2=folded_operator_tensor(C_array,O2,D,chi)
            fold_c=foldedchannel(layer_to_moves[Lx-1],D,chi)
            start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[1,2,3,-1],[2],[3]])
            fold_c=foldedchannel(C_array,D,chi)
            for i in range(y-1):
                start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[2,1,3,-1],[2],[3]])
            return ncon([start_c,fold_o2,fold_in,fold_in,fold_in],[[1],[2,1,3,4],[2],[3],[4]])
        '''
    elif y2>=y1 and x2>x1:
        pass
        
    else:
        print("Not considered yet\n")
        exit(1)

def expectation_value(fold_c,fold_o,y,chi):
    fold_in=folded_identity(chi)
    if y==0:
        return ncon([fold_o,fold_in,fold_in,fold_in,fold_in],[[1,2,3,4],[1],[2],[3],[4]])
    start_c=ncon([fold_c,fold_in,fold_in,fold_in],[[1,2,3,-1],[1],[2],[3]])
    for i in range(y-1):
        start_c=ncon([start_c,fold_c,fold_in,fold_in],[[1],[2,1,3,-1],[2],[3]])
    return ncon([start_c,fold_o,fold_in,fold_in,fold_in],[[1],[2,1,3,4],[2],[3],[4]])

def expectation_value_from_op(O,y,C_array,D,chi):
    fold_o=folded_operator_tensor(C_array,O,D,chi)
    fold_c=foldedchannel(C_array,D,chi)
    return expectation_value(fold_c,fold_o,y,chi)

def test_expectation_value(D,chi,y):
    op=np.random.random([D,D])+1j*np.random.random([D,D])
    from decomposition import random_dual_PEPS
    C_array=random_dual_PEPS(D=D,chi=chi)
    temp1=expectation_value_from_op(op,y,C_array,D,chi)
    temp2=general_MPS_direct_contraction_from_op(np.eye(D),op,0,0,1,y,C_array,D,chi)
    print(temp1)
    print(temp2)
    print(np.abs(temp1-temp2)/np.abs(temp1))


def check_MPS_PEPS_coincident(D,chi,L):
    from decomposition import random_dual_PEPS
    #C_array=random_dual_PEPS(D,chi)
    U=dual_sequential_gates(D)
    C_array=construct_equal_dual_ISOTNS_from_gate(U,chi)
    #U=controlled_dual_unitary_gates(D)
    #C_array=construct_equal_dual_ISOTNS_from_gate(U,chi)
    #U=sequential_2qubitgates()
    #C_array=construct_equal_dual_ISOTNS_from_gate(U,chi)
    #U=translation_gates(D)
    #C_array=construct_equal_dual_ISOTNS_from_gate(U,chi)
    O1=np.random.random([D,D])+1j*np.random.random([D,D])
    O2=np.random.random([D,D])+1j*np.random.random([D,D])
    O1=O1+O1.conj().T
    O2=O2+O2.conj().T
    O1=O1-np.trace(O1)*np.eye(D)/D
    O2=O2-np.trace(O2)*np.eye(D)/D
    #O1=np.eye(D)
    y1=np.random.randint(0,L)
    y2=np.random.randint(y1,L)
    x1=np.random.randint(0,L) 
    while (x1==L-1):
        y1=np.random.randint(0,L)
        y2=np.random.randint(y1,L)
        x1=np.random.randint(0,L) 
    x2=np.random.randint(x1+1,L)
    y1=0
    independent_multiple=expectation_value_from_op(O1,y1,C_array,D,chi)*expectation_value_from_op(O2,y2,C_array,D,chi)
    temp1=PEPS_direct_from_op(O1,O2,x1,y1,x2,y2,C_array,D,chi,L,L)
    temp2=general_MPS_direct_contraction_from_op(O1,O2,x1,y1,x2,y2,C_array,D,chi)
    temp_normalize=general_MPS_direct_contraction_from_op(np.eye(D),np.eye(D),x1,y1,x2,y2,C_array,D,chi)
    print("\n\n\n")
    print("({},{})   ({},{})".format(x1,y1,x2,y2))
    print("The exact result, unnormalized result, normalization factor, normalized result, independent multiply and normalized error are:")
    print(temp1-independent_multiple)
    print(temp2-independent_multiple)
    print(temp_normalize)
    print(temp2/temp_normalize-independent_multiple)
    print(independent_multiple)
    print(abs(temp1-temp2/temp_normalize)/abs(temp1))

def check_from_folded_channel(D,chi,L):
    from decomposition import random_dual_PEPS_fold
    fold_c=random_dual_PEPS_fold(D,chi)
    fold_o1=np.random.random([chi**2]*4)+1j*np.random.random([chi**2]*4)
    fold_o2=np.random.random([chi**2]*4)+1j*np.random.random([chi**2]*4)
    O1=O1+O1.conj().T
    O2=O2+O2.conj().T
    O1=O1-np.trace(O1)*np.eye(D)/D
    O2=O2-np.trace(O2)*np.eye(D)/D
    x1=np.random.randint(0,L)
    y1=0
    x2=np.random.randint(0,L)
    y2=np.random.randint(0,L)
    while x1==x2 and y1==y2:
        y2=np.random.randint(0,L)
    temp1=MPS_direct_contraction(fold_o1,fold_o2,x2-x1,y2-y1,fold_c,chi)
    temp2=PEPS_direct_contraction(fold_o1,fold_o2,x1,y1,x2,y2,fold_c,chi,L,L)
    print("\n\n\n")
    print("({},{})   ({},{})".format(x1,y1,x2,y2))
    print(temp1)
    print(temp2)
          


if __name__=='__main__':
    check_MPS_PEPS_coincident(2,2,6)