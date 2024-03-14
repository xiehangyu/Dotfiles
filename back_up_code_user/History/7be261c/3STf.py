import numpy as np

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
    CM=np.zeros([BD**2*2,PD*BD**2*2])
    condi_in=0
    for i in range(BD):
        for j in range(i,BD):
            vari_in=0
            if i!=j:
                for n in range(BD):
                    for m in range(PD):
                        for k in range(BD):
                            CM[condi_in,vari_in]=np.imag(R_ij_C_nmk_re(c_array,i,j,n,m,k))
                            CM[condi_in+1,vari_in]=np.imag(L_ij_C_nmk_re(c_array,i,j,n,m,k))
                            vari_in+=1
                            CM[condi_in,vari_in]=np.imag(R_ij_C_nmk_im(c_array,i,j,n,m,k))

                        