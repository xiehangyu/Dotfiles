import numpy as np
import matplotlib.pyplot as plt
from contraction import *
from test_topolophase import *


def folded_W_doubled(W):
    W=np.reshape(W,[2,2,2,2])
    W_double=np.zeros([2]*8)
    W_double=ncon([W,W],[[-1,-3,-5,-7],[-2,-4,-6,-8]])
    return W_double.reshape([4,4,4,4])


def reduced_W_expect(W,folded_o,y):
    W1=folded_W_doubled(W)
    W2=foldedchannel(construct_C_array_from_W(W),16,2)
    fold_in=folded_identity(2)
    if y==0:
        start=ncon([W1,fold_in,fold_in,fold_in],[[[1,2,-1,3],[1],[2],[3]]])
        start=ncon([start,folded_o,W1,fold_in,fold_in,fold_in,fold_in],[[1],[1,2,3,-1],[3,4,5,6],[2],[4],[5],[6]])
        return ncon([start,W1,fold_in,fold_in,fold_in],[[1],[2,1,3,4],[2],[3],[4]])
    else:
        start=ncon([W1,folded_o,W1,W1]+[fold_in]*7,[[1,-1,3,2],[3,-2,5,4],[5,-3,7,6],[8,4,10,9],[1],[2],[6],[7],[8],[10],[9]])
        start=ncon([start,W2,W1,W2,fold_in,fold_in],[[1,2,3],[4,-1,5,1],[5,-2,6,2],[6,-3,7,3],[4],[7]])
        if y==1:
            return ncon([start,fold_in,fold_in,fold_in],[[1,2,3],[1],[2],[3]])
        for i in range(1,y-1):
            start=ncon([start,W2,W2,W2,fold_in,fold_in],[[1,2,3],[4,-1,5,1],[5,-2,6,2],[6,-3,7,3],[4],[7]])
        return ncon([start,fold_in,fold_in,fold_in],[[1,2,3],[1],[2],[3]])
    
def reduced_correlator(fold_O1,fold_O2,x,y,W):
    #Here we assume that both x,y>=3
    W1=folded_W_doubled(W)
    W2=foldedchannel(construct_C_array_from_W(W),16,2)
    fold_in=folded_identity(2)
    left_start=ncon([W1,fold_in,fold_in,fold_in],[[1,2,-1,3],[1],[2],[3]])
    left_start=ncon([left_start,fold_O1,W1,fold_in,fold_in,fold_in],[[1],[1,2,-1,4],[3,4,-2,5],[2],[3],[5]])
    left_start=ncon([left_start,W1,W2,fold_in,fold_in],[[2,1],[2,3,-1,4],[1,4,-2,5],[3],[5]])
    up_start=ncon([W1,fold_O2,W1,W1]+[fold_in]*7,[[1,-1,3,2],[3,-2,5,4],[5,-3,7,6],[8,4,10,9],[1],[2],[6],[7],[8],[10],[9]])
    up_start=ncon([up_start,W2,W1,W2,fold_in,fold_in],[[1,2,3],[4,-1,5,1],[5,-2,6,2],[6,-3,7,3],[4],[7]])
    for i in range(x-3):
        left_start=ncon([left_start,W2,W2,fold_in,fold_in],[[1,2],[1,3,-1,4],[2,4,-2,5],[3],[5]])
    for i in range(y-3):
        up_start=ncon([up_start,W2,W2,W2,fold_in,fold_in],[[1,2,3],[4,-1,5,1],[5,-2,6,2],[6,-3,7,3],[4],[7]])
    rd_start=ncon([W2,W2,W2,fold_in,fold_in,fold_in,fold_in],[[-1,1,2,-2],[2,3,4,-3],[4,5,6,-4],[1],[3],[5],[6]])
    rd_start=ncon([rd_start,W2,W2,W2,fold_in],[[-1,1,2,3],[-2,1,4,-3],[4,2,5,-4],[5,3,6,-5],[6]])
    return ncon([left_start,up_start,rd_start],[[1,2],[3,4,5],[1,2,3,4,5]])


def construct_C_array_from_W(W):
    new_W=W.reshape([2,2,2,2])
    C_array=np.zeros([2]*8)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    C_array[i,j,k,l,i,j,k,l]=new_W[i,j,k,l]
    return C_array.reshape([16,2,2,2,2])


def reduced_calculate_correlators(W,positions,op1=None,op2=None):
    #We assume that x1=y1=0
    C_array=construct_C_array_from_W(W)
    if op1==None:
        op1=np.random.random([16,16])+1j*np.random.random([16,16])
        op2=np.random.random([16,16])+1j*np.random.random([16,16])
        op1=op1+np.conj(op1).T
        op2=op2+np.conj(op2).T
        op1=op1-np.trace(op1)*np.eye(16)/16
        op2=op2-np.trace(op2)*np.eye(16)/16
        op1=op1*1000
        op2=op2*1000
    fold_O1=folded_operator_tensor(C_array,op1,16,2)
    fold_O2=folded_operator_tensor(C_array,op2,16,2)
    expect1=reduced_W_expect(W,fold_O1,0)
    print("The expectation value 1 is {}".format(expect1))
    correlators=[]
    for x,y in positions:
        print("({},{})".format(x,y))
        expect2=reduced_W_expect(W,fold_O2,y)
        print("The expectation value 2 is {}".format(expect2))
        cor=MPS_direct_contraction(fold_O1,fold_O2,x,y,fold_c,2)
        print("The correlator is {}".format(cor))
        correlators.append(cor-expect1*expect2)
    return correlators

def calculate_correlators(W,positions,op1=None,op2=None):
    #We assume that x1=y1=0
    C_array=construct_C_array_from_W(W)
    fold_c=foldedchannel(C_array,16,2)
    if op1==None:
        op1=np.random.random([16,16])+1j*np.random.random([16,16])
        op2=np.random.random([16,16])+1j*np.random.random([16,16])
        op1=op1+np.conj(op1).T
        op2=op2+np.conj(op2).T
        op1=op1-np.trace(op1)*np.eye(16)/16
        op2=op2-np.trace(op2)*np.eye(16)/16
    fold_O1=folded_operator_tensor(C_array,op1,16,2)
    fold_O2=folded_operator_tensor(C_array,op2,16,2)
    expect1=expectation_value(fold_c,fold_O1,0,2)
    print("The expectation value 1 is {}".format(expect1))
    correlators=[]
    for x,y in positions:
        print("({},{})".format(x,y))
        expect2=expectation_value(fold_c,fold_O2,y,2)
        print("The expectation value 2 is {}".format(expect2))
        cor=MPS_direct_contraction(fold_O1,fold_O2,x,y,fold_c,2)
        print("The correlator is {}".format(cor))
        correlators.append(cor-expect1*expect2)
    return correlators


def drawing():
    #The horizontal second largest eigenvalue should be |alpha+beta-1|,eigenstate is (1,-1)
    #The vertical second largest eigenvalue should be |alpha+beta-1|,eigenstate is (1,-1)
    #The vertical largest eigenvalue should be 1, the eigenvector is (1-beta,1-alpha)
    W=create_W_matrix3(0.3,0.4)
    positions=[]
    ilist=[]
    for i in range(1,20):
        positions.append((0,i))
        ilist.append(i)
    correlators=calculate_correlators(W,positions)
    print(correlators)
    plt.plot(ilist,correlators)
    plt.show()
    correlators=np.log(np.abs(np.array(correlators)))
    plt.plot(ilist,correlators)
    print("The decay rate is {}".format(np.exp((correlators[11]-correlators[3])/(ilist[11]-ilist[3]))))
    plt.show()
    positions=[]
    ilist=[]
    for i in range(1,20):
        positions.append((i,0))
        ilist.append(i)
    correlators=calculate_correlators(W,positions)
    print(correlators)
    plt.plot(ilist,correlators)
    plt.show()
    correlators=np.log(np.abs(np.array(correlators)))
    plt.plot(ilist,correlators)
    print("The decay rate is {}".format(np.exp((correlators[15]-correlators[3])/(ilist[15]-ilist[3]))))
    plt.show()
    positions=[]
    ilist=[]
    for i in range(1,20):
        positions.append((i,i))
        ilist.append(i)
    correlators=calculate_correlators(W,positions)
    print(correlators)
    plt.plot(ilist,correlators)
    plt.show()
    correlators=np.log(np.abs(np.array(correlators)))
    plt.plot(ilist,correlators)
    print("The decay rate is {}".format(np.exp((correlators[15]-correlators[3])/(ilist[15]-ilist[3]))))
    plt.show()


if __name__=='__main__':
    drawing()