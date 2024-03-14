import numpy as np
import sys
from scipy.stats import unitary_group
from scipy.linalg import expm

def R_ij_C_nmk(c_array,i,j,n,m,k):
    if j!=n:
        return 0
    else:
        return np.conjugate(c_array[i,m,k])

def R_ij_Cconj_nmk(c_array,i,j,n,m,k):
    if i!=n:
        return 0
    else:
        return c_array[j,m,k]

def L_ij_C_nmk(c_array,i,j,n,m,k):
    if j!=k:
        return 0
    else:
        return np.conjugate(c_array[n,m,i])


def L_ij_Cconj_nmk(c_array,i,j,n,m,k):
    if i!=k:
        return 0
    else:
        return c_array[n,m,i]

def R_ij_C_nmk_re(c_array,i,j,n,m,k):
    return R_ij_C_nmk(c_array,i,j,n,m,k)+R_ij_Cconj_nmk(c_array,i,j,n,m,k)

def R_ij_C_nmk_im(c_array,i,j,n,m,k):
    return 1j*(R_ij_C_nmk(c_array,i,j,n,m,k)-R_ij_Cconj_nmk(c_array,i,j,n,m,k))

def L_ij_C_nmk_re(c_array,i,j,n,m,k):
    return L_ij_C_nmk(c_array,i,j,n,m,k)+L_ij_Cconj_nmk(c_array,i,j,n,m,k)

def L_ij_C_nmk_im(c_array,i,j,n,m,k):
    return 1j*(L_ij_C_nmk(c_array,i,j,n,m,k)-L_ij_Cconj_nmk(c_array,i,j,n,m,k))


def isoMPS_conditionmatrix(c_array,PD,BD):
    CM=np.zeros([BD**2,PD*BD**2*2],dtype=np.float64)
    condi_in=0
    for i in range(BD):
        for j in range(i,BD):
            if i!=j:
                vari_in=0
                for n in range(BD):
                    for m in range(PD):
                        for k in range(BD):
                            CM[condi_in,vari_in]=np.imag(L_ij_C_nmk_re(c_array,i,j,n,m,k))
                            vari_in+=1
                            CM[condi_in,vari_in]=np.imag(L_ij_C_nmk_im(c_array,i,j,n,m,k))
                            vari_in+=1
                condi_in+=1
            vari_in=0
            for n in range(BD):
                for m in range(PD):
                    for k in range(BD):
                        CM[condi_in,vari_in]=np.real(L_ij_C_nmk_re(c_array,i,j,n,m,k))
                        vari_in+=1
                        CM[condi_in,vari_in]=np.real(L_ij_C_nmk_im(c_array,i,j,n,m,k))
                        vari_in+=1
            condi_in+=1
    print("The condition matrix is ({},{})".format(condi_in,vari_in))
    return CM

def MPS_conditionmatrix(c_array,PD,BD):
    CM=np.zeros([BD**2*2,PD*BD**2*2],dtype=np.float64)
    condi_in=0
    for i in range(BD):
        for j in range(i,BD):
            if i!=j:
                vari_in=0
                for n in range(BD):
                    for m in range(PD):
                        for k in range(BD):
                            CM[condi_in,vari_in]=np.imag(R_ij_C_nmk_re(c_array,i,j,n,m,k))
                            CM[condi_in+1,vari_in]=np.imag(L_ij_C_nmk_re(c_array,i,j,n,m,k))
                            vari_in+=1
                            CM[condi_in,vari_in]=np.imag(R_ij_C_nmk_im(c_array,i,j,n,m,k))
                            CM[condi_in+1,vari_in]=np.imag(L_ij_C_nmk_im(c_array,i,j,n,m,k))
                            vari_in+=1
                condi_in+=2
            vari_in=0
            for n in range(BD):
                for m in range(PD):
                    for k in range(BD):
                        CM[condi_in,vari_in]=np.real(R_ij_C_nmk_re(c_array,i,j,n,m,k))
                        CM[condi_in+1,vari_in]=np.real(L_ij_C_nmk_re(c_array,i,j,n,m,k))
                        vari_in+=1
                        CM[condi_in,vari_in]=np.real(R_ij_C_nmk_im(c_array,i,j,n,m,k))
                        CM[condi_in+1,vari_in]=np.real(L_ij_C_nmk_im(c_array,i,j,n,m,k))
                        vari_in+=1
            condi_in+=2
    print("The condition matrix is ({},{})".format(condi_in,vari_in))
    return CM


def check_isoTNS_dimension(BD):
    total_param=BD**3*2
    inte_matrix=unitary_group.rvs(BD**2)
    c_array=construct_c_array_equalPB(inte_matrix,BD)
    CM=isoMPS_conditionmatrix(c_array,BD,BD)
    rank=np.linalg.matrix_rank(CM)
    print("The tangent dimension is {}".format(total_param-rank))

def tangent_dimension_local(c_array,PD,BD):
    total_param=PD*BD**2*2
    CM=MPS_conditionmatrix(c_array,PD,BD)
    rank=np.linalg.matrix_rank(CM)
    print("The eigenvalues are {}".format(np.linalg.eigvals(np.dot(CM,CM.conjugate().transpose()))))
    return total_param-rank
    
def construct_c_array_equalPB(intm_mattrix,BD):
    c_array=np.zeros([BD,BD,BD],dtype=np.complex128)
    for n in range(BD):
        for m in range(BD):
            for k in range(BD):
                c_array[n,m,k]=intm_mattrix[m*BD+k,n*BD]
    return c_array


def check_c_array(c_array,PD,BD):
    temp1=c_array.reshape(BD*PD,BD)
    identity=np.eye(BD)
    print("The left canonicality is {}".format(np.linalg.norm(np.dot(temp1.conjugate().transpose(),temp1)-identity)))
    temp2=c_array.reshape(BD,BD*PD)
    print("The right canonicality is {}".format(np.linalg.norm(np.dot(temp2,temp2.conjugate().transpose())-identity)))

def general_SWAP_operator(D):
    SWAP=np.zeros([D**2,D**2],dtype=np.complex128)
    for i in range(D):
        for j in range(D):
            SWAP[i*D+j,j*D+i]=1
    return SWAP


def random_general_dual_unitary(D):
    u1=unitary_group.rvs(D)
    u2=unitary_group.rvs(D)
    u3=unitary_group.rvs(D)
    u4=unitary_group.rvs(D)
    j1=np.random.random()
    SWAP=general_SWAP_operator(D)
    Z_general=np.zeros([D,D],dtype=np.complex128)
    for i in range(D):
        Z_general[i,i]=D/2-i
    ZZ_general=np.kron(Z_general,Z_general)
    inte_matrix=np.kron(u1,u2).dot(SWAP).dot(expm(2j*j1*np.pi*ZZ_general)).dot(np.kron(u3,u4))
    return inte_matrix
    


def first_trial():
    SWAP=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],dtype=np.complex128)
    Z=np.array([[1,0],[0,-1]],dtype=np.complex128)
    ZZ=np.kron(Z,Z)
    X=np.array([[0,1],[1,0]],dtype=np.complex128)
    XX=np.kron(X,X)
    Y=np.array([[0,-1j],[1j,0]],dtype=np.complex128)
    YY=np.kron(Y,Y)
    j1=np.random.random()
    j2=np.random.random()
    u1=unitary_group.rvs(2)
    u2=unitary_group.rvs(2)
    u3=unitary_group.rvs(2)
    inte_matrix=np.kron(u1,u2).dot(expm(1j*j1*np.pi*XX)).dot(expm(1j*j2*np.pi*ZZ)).dot(expm(1j*np.pi/4*YY)).dot(np.kron(u3,np.eye(2)))
    c_array=construct_c_array_equalPB(inte_matrix,2)
    check_c_array(c_array,2,2)
    print("The tangent dimension is {}".format(tangent_dimension_local(c_array,2,2)))


def second_trial():
    SWAP=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],dtype=np.complex128)
    Z=np.array([[1,0],[0,-1]],dtype=np.complex128)
    ZZ=np.kron(Z,Z)
    j1=np.random.random()
    j1=1.0/4.0
    u1=unitary_group.rvs(2)
    u2=unitary_group.rvs(2)
    u3=unitary_group.rvs(2)
    u4=unitary_group.rvs(2)
    inte_matrix=np.kron(u1,u2).dot(SWAP).dot(expm(1j*j1*np.pi*ZZ)).dot(np.kron(u3,u4))
    c_array=construct_c_array_equalPB(inte_matrix,2)
    check_c_array(c_array,2,2)
    print("The tangent dimension is {}".format(tangent_dimension_local(c_array,2,2)))

def third_trial(D):
    inte_matrix=random_general_dual_unitary(D)
    c_array=construct_c_array_equalPB(inte_matrix,D)
    check_c_array(c_array,D,D)
    print("The tangent dimension is {}".format(tangent_dimension_local(c_array,D,D)))


if __name__=="__main__":
    #second_trial()
    third_trial(int(sys.argv[1]))
     
                        