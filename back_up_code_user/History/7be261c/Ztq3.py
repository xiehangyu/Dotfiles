import numpy as np
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


def tangent_dimension_local(c_array,PD,BD):
    total_param=PD*BD**2*2
    CM=MPS_conditionmatrix(c_array,PD,BD)
    rank=np.linalg.matrix_rank(CM)
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

def first_trial():
    SWAP=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],dtype=np.complex128)
    Z=np.array([[1,0],[0,-1]],dtype=np.complex128)
    ZZ=np.kron(Z,Z)
    inte_matrix=np.dot(SWAP,expm(1j*ZZ))
    c_array=construct_c_array_equalPB(inte_matrix,2)
    check_c_array(c_array,2,2)
    print("The tangent dimension is {}".format(tangent_dimension_local(c_array,2,2)))

if __name__=="__main__":
    first_trial()
     
                        