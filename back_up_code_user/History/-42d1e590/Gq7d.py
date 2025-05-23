from ncon import ncon
import math
import numpy as np
import numpy.linalg as LA
from contraction import *
from decomposition import *
import scipy
import matplotlib.pyplot as plt
from calculate_dual_ISOTNS import *

def contract_identity_n_depth(fold_c,Lx,Ly,chi):
    fold_in=folded_identity(chi)
    start_c=ncon([fold_c,fold_in,fold_in],[[1,2,-1,-2],[1],[2]])
    for i in range(1,Lx):
        start_c=ncon([start_c,fold_c,fold_in],[[1]+[-j for j in range(3,i+3)],[1,2,-1,-2],[2]])
    start_c=ncon([start_c,fold_in],[[1]+[-j for j in range(1,Lx+1)],[1]])
    for i in range(1,Ly):
        start_c=ncon([start_c,fold_c,fold_in],[[-j for j in range(1,Lx)]+[2],[1,2,-Lx,-Lx-1],[1]])
        for j in range(1,Lx):
            start_c=ncon([start_c,fold_c],[[-k for k in range(1,Lx-j)]+[2,1]+[-k for k in range(Lx-j+2,Lx+2)],[1,2,-Lx+j,-Lx+j-1]])
        start_c=ncon([start_c,fold_in],[[1]+[-j for j in range(1,Lx+1)],[1]])
    return start_c

def testing_contract_identity(Lx,Ly,D=2,chi=2):
    fold_c=random_dual_PEPS_fold(D,chi)
    final_c=contract_identity_n_depth(fold_c,Lx,Ly,chi)
    fold_in=folded_identity(chi)
    return ncon([final_c]+[fold_in for i in range(Lx)],[[j for j in range(1,Lx+1)]]+[[j] for j in range(1,Lx+1)])

def distance_for_multiply(fold_c,Lx,Ly,chi):
    fold_in=folded_identity(chi)
    reference_state=ncon([fold_in]*Lx,[[-j] for j in range(1,Lx+1)])
    final_c=contract_identity_n_depth(fold_c,Lx,Ly,chi)
    return LA.norm(final_c/LA.norm(final_c)-reference_state)

def is_normal(A):
    return LA.norm(A@A.conj().T-A.conj().T@A)/(LA.norm(A)*LA.norm(A.conj().T))




def rounded_N_matrix(N,Lx,chi):
    states=[]
    for i in range(N):
        fold_c=random_dual_PEPS_fold(chi,chi)
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
        states.append(eigenstate)
    


def test_eigenvalues_entropy(fold_c,Lx,chi):
    htm=horizontal_transform_matrix(fold_c,Lx,chi)
    bd=chi**2
    htm=htm.transpose([j for j in range(Lx)]+[2*Lx-j-1 for j in range(Lx)])
    htm=htm.reshape([bd**Lx]*2).transpose()
    normal_difference=is_normal(htm)
    eigenvalues,eigenstates=scipy.linalg.eig(htm)
    idx=eigenvalues.argsort()[::-1]
    eigenvalues=eigenvalues[idx]
    eigenstates=eigenstates[:,idx]
    print(eigenvalues[:5])
    eigenstate=eigenstates[:,0]
    reference_eigen_vector=kron_product([folded_identity(chi) for i in range(Lx)])
    distance_temp=LA.norm(eigenstate-reference_eigen_vector)
    print("The difference between the eigenstate and the reference state is ",distance_temp)
    print("Another way to calculate distance is {}".format(distance_for_multiply(fold_c,Lx,12,chi)))
    dx=bd**(Lx//2)
    eigenstateM=eigenstate.reshape([dx,bd**(Lx-Lx//2)])
    eigenstateM=eigenstateM.dot(eigenstateM.conj().T)
    eigenvaluesM=LA.eigvalsh(eigenstateM)
    entropy=0
    print("The trace is :")
    print(np.sum(eigenvaluesM))
    for i in range(len(eigenvaluesM)):
        if eigenvaluesM[i]<=0 and eigenvaluesM[i]>=-1e-10:
            continue
        elif eigenvaluesM[i]>0:
            entropy+=-eigenvaluesM[i]*np.log(eigenvaluesM[i])
        else:
            print("eigenvalue is negative")
            print(eigenvaluesM)
            return False
    print("The entropy is ",entropy)
    return entropy,distance_temp,normal_difference

if __name__=='__main__':
    #fold_c=random_dual_PEPS_fold(16,2)
    gate=controlled_dual_unitary_gates(2)
    C_array=construct_equal_dual_ISOTNS_from_gate(gate,2)
    fold_c=foldedchannel(C_array,2,2)
    Ls=[]
    entropys=[]
    distances=[]
    normal_difs=[]
    for i in [2,3,4,5,6]:
        Ls.append(i)
        entro,dista,norm_dif=test_eigenvalues_entropy(fold_c,i,2)
        entropys.append(entro)
        distances.append(dista)
        normal_difs.append(norm_dif)
    plt.plot(Ls,entropys)
    plt.show()
    plt.plot(Ls,distances)
    plt.xlabel("Lx")
    plt.ylabel("distance")
    plt.show()
    plt.plot(Ls,normal_difs)
    plt.xlabel("Lx")
    plt.ylabel("regularized norm of commutator")
    plt.show()
    