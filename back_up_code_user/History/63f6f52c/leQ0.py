import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.linalg import expm
from finding import *
from scipy.stats import unitary_group
def Clifford_group(D,p,q):
    translation_operator=np.zeros([D,D],dtype=complex)
    for i in range(D):
        translation_operator[(i+p)%D,i]=1
    phase_operator=np.zeros([D,D],dtype=complex)
    omega=np.exp(2*np.pi*1j/D)
    for i in range(D):
        phase_operator[i,i]=omega**(i*q)
    return translation_operator.dot(phase_operator)

def matrixtovector(M):
    vec=np.zeros(np.shape(M)[0]**2,dtype=complex)
    for i in range(np.shape(M)[0]):
        for j in range(np.shape(M)[1]):
            vec[i*np.shape(M)[0]+j]=M[i,j]
    return np.array([list(vec)])

def generate_matrix(D,alpha,beta,gamma):
    return_matrix=np.zeros([D**2,D**2],dtype=complex)
    omega=np.exp(2*np.pi*1j/D)
    for p in range(D):
        for q in range(D):
            phase=omega**(alpha*p*p+beta*p*q+gamma*q*q)
            temp_matrix=Clifford_group(D,p,q)
            temp_vec=matrixtovector(temp_matrix)
            return_matrix += phase*temp_vec.conj().transpose().dot(temp_vec)
    return return_matrix/D

def epsilonL(U,a):
    U=np.copy(U.conj().transpose())
    dim=np.shape(a)[0]
    temp_matrix=np.kron(np.eye(dim),a)
    temp_matrix=U.conj().transpose().dot(temp_matrix.dot(U))
    return_matrix=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim))
        return_matrix += project.dot(temp_matrix.dot(project.conj().transpose()))
    return return_matrix

def epsilonR(U,a):
    U=np.copy(U.conj().transpose())
    dim=np.shape(a)[0]
    temp_matrix=np.kron(a,np.eye(dim))
    temp_matrix=U.conj().transpose().dot(temp_matrix.dot(U))
    return_matrix=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim),new_vector)
        return_matrix += project.dot(temp_matrix.dot(project.conj().transpose()))
    return return_matrix

def ML(U,a):
    U=np.copy(U.conj().transpose())
    dim=np.shape(a)[0]
    temp_matrix=np.kron(a,np.eye(dim))
    temp_matrix=U.conj().transpose().dot(temp_matrix.dot(U))
    return_matrix=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim))
        return_matrix += project.dot(temp_matrix.dot(project.conj().transpose()))
    return return_matrix

def MR(U,a):
    U=np.copy(U.conj().transpose())
    dim=np.shape(a)[0]
    temp_matrix=np.kron(np.eye(dim),a)
    temp_matrix=U.conj().transpose().dot(temp_matrix.dot(U))
    return_matrix=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim),new_vector)
        return_matrix += project.dot(temp_matrix.dot(project.conj().transpose()))
    return return_matrix

def interaction(j1,j2,j3):
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    XX=np.kron(X,X)
    YY=np.kron(Y,Y)
    ZZ=np.kron(Z,Z)
    return expm(-1j*(j1*XX+j2*YY+j3*ZZ))

def singleoperator(r,theta,phi):
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    H=r*(np.cos(theta)*Z+np.sin(theta)*np.cos(phi)*X+np.sin(theta)*np.sin(phi)*Y)
    return expm(1j*H)


def general_controlZoperator(D):
    phase_operator=np.zeros((D,D),dtype=complex)
    for i in range(D):
        phase_operator[i,i]=np.exp(2*np.pi*1j*i/D)
    return_operator=np.zeros((D**2,D**2),dtype=complex)
    for i in range(D):
        project_op=np.zeros((D,D),dtype=complex)
        project_op[i,i]=1
        return_operator += np.kron(project_op,np.linalg.matrix_power(phase_operator,i))
    return return_operator



def general_control_operator(D):
    shift_matrix=np.zeros((D,D),dtype=complex)
    for i in range(D):
        shift_matrix[(i+1)%D,i]=1
    return_operator=np.zeros((D**2,D**2),dtype=complex)
    for i in range(D):
        project_op=np.zeros((D,D),dtype=complex)
        project_op[i,i]=1
        return_operator += np.kron(project_op,np.linalg.matrix_power(shift_matrix,i))
    return return_operator

def channeltwosite(U,a):
    U=np.copy(U.conj().transpose())
    dim=int(np.sqrt(np.shape(a)[0]))
    temp_matrix=np.kron(np.eye(dim),a)
    temp_matrix=np.kron(temp_matrix,np.eye(dim))
    temp_matrix=np.kron(U.conj().transpose(),U.conj().transpose()).dot(temp_matrix.dot(np.kron(U,U)))
    return_matrix=np.zeros((dim**3,dim**3),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim**3))
        return_matrix += project.dot(temp_matrix.dot(project.conj().transpose()))
    temp_matrix=np.copy(return_matrix)
    return_matrix=np.zeros((dim**2,dim**2),dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim**2),new_vector)
        return_matrix += project.dot(temp_matrix.dot(project.conj().transpose()))
    return return_matrix

def twositeoperator1leftstep(t,pa,pb,U):
    origin_t=t
    dim=int(np.sqrt(np.shape(pa)[0]))
    pb=U.dot(pb.dot(U.conj().transpose()))
    newpb=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim),new_vector)
        newpb += project.dot(pb.dot(project.conj().transpose()))
    t=t-1
    pb=newpb
    while t>1:
        pb=epsilonL(U,pb)
        pb=epsilonR(U,pb)
        t=t-2
    if t==0:
        pb=np.kron(np.eye(dim),pb)
        return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))
    if t==1:
        origin_t=origin_t-1
        pb=np.kron(np.eye(dim),pb)
        pb=U.dot(pb.dot(U.conj().transpose()))
        return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def twositeoperator1rightstep(t,pa,pb,U):
    dim=int(np.sqrt(np.shape(pa)[0]))
    origin_t=t
    pb=U.dot(pb.dot(U.conj().transpose()))
    newpb=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim))
        newpb += project.dot(pb.dot(project.conj().transpose()))
    t=t-1
    pb=newpb
    while t>1:
        pb=epsilonR(U,pb)
        pb=epsilonL(U,pb)
        t=t-2
    if t==0:
        pb=np.kron(pb,np.eye(dim))
        return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))
    if t==1:
        origin_t=origin_t-1
        pb=np.kron(pb,np.eye(dim))
        pb=U.dot(pb.dot(U.conj().transpose()))
        return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def twositeoperatorleftlightcone(t,pa,pb,U):
    origin_t=t
    dim=int(np.sqrt(np.shape(pa)[0]))
    pb=U.dot(pb.dot(U.conj().transpose()))
    newpb=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim),new_vector)
        newpb += project.dot(pb.dot(project.conj().transpose()))
    t=t-1
    pb=newpb
    while t>1:
        pb=ML(U,pb)
        t=t-1
    pb=np.kron(np.eye(dim),pb)
    pb=U.dot(pb.dot(U.conj().T))
    origin_t=origin_t-1
    return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def twositeoperatorrightlightcone(t,pa,pb,U):
    origin_t=t
    dim=int(np.sqrt(np.shape(pa)[0]))
    pb=U.dot(pb.dot(U.conj().transpose()))
    newpb=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim))
        newpb += project.dot(pb.dot(project.conj().transpose()))
    t=t-1
    pb=newpb
    while t>1:
        pb=MR(U,pb)
        t=t-1
    pb = np.kron(pb,np.eye(dim))
    pb = U.dot(pb.dot(U.conj().T))
    origin_t=origin_t-1
    return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def twositeoperatorleftlightcone1step(t,pa,pb,U):
    dim=int(np.sqrt(np.shape(pa)[0]))
    newpa=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim))
        newpa += project.dot(pa.dot(project.conj().transpose()))
    newpa=newpa/dim
    newpa=np.kron(newpa,np.eye(dim))
    return twositeoperatorleftlightcone(t,newpa,pb,U)

def twositeoperatorrightlightcone1step(t,pa,pb,U):
    dim=int(np.sqrt(np.shape(pa)[0]))
    newpa=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim),new_vector)
        newpa += project.dot(pa.dot(project.conj().transpose()))
    newpa=newpa/dim
    newpa=np.kron(np.eye(dim),newpa)
    return twositeoperatorrightlightcone(t,newpa,pb,U)
        

def twositeoperator(t,pa,pb,U1,U2):
    U1=np.copy(U1.conj().transpose())
    dim=np.shape(pa)[0]
    origin_t=0
    while t>0:
        pb=U1.conj().transpose().dot(pb.dot(U1))
        t=t-1
        origin_t+=1
        if t>0:
            pb=channeltwosite(U2,pb)
            t=t-1
    return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def time_correlation_function(t,L,pa,pb,U1,U2):
    dim=np.shape(pa)[0]
    origin_t=t
    if L==0:
        while t>0:
            pb=epsilonR(U1,pb)
            t=t-1
            if t>0:
                pb=epsilonL(U2,pb)
                t=t-1
    else:
        while t>0:
            pb=epsilonL(U2,pb)
            t=t-1
            if t>0:
                pb=epsilonR(U1,pb)
                t=t-1
    return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def support2correlationfunction(t,j,i,pa,pb,U):
    if i%2!=0:
        print("i must be even")
        return
    if t==0:
        return np.trace(pa.dot(pb))
    if j==i:
        return twositeoperator(t,pa,pb,U,U)
    if j==i-1:
        return twositeoperator1leftstep(t,pa,pb,U)
    if j==i+1:
        return twositeoperator1rightstep(t,pa,pb,U)
    if j==i-t+1:
        return twositeoperatorleftlightcone(t,pa,pb,U)
    if j==i-t:
        return twositeoperatorleftlightcone1step(t,pa,pb,U)
    if j==i+t-1:
        return twositeoperatorrightlightcone(t,pa,pb,U)
    if j==i+t:
        return twositeoperatorrightlightcone1step(t,pa,pb,U)
    return 0

def lightconechannelthreesite(U,a):
    D=int(np.sqrt(np.shape(a)[0]))
    anew=np.kron(a,np.eye(D))
    U1=np.kron(np.eye(D),U)
    U2=np.kron(U,np.eye(D))
    anew=U1.dot(anew.dot(U1.conj().transpose()))
    anew=U2.dot(anew.dot(U2.conj().transpose()))
    return_matrix=np.zeros((D**2,D**2),dtype=complex)
    zero_vector=np.zeros(D,dtype=complex)
    for i in range(D):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(D**2))
        return_matrix += project.dot(anew.dot(project.conj().transpose()))
    return return_matrix

def cf3sitesrightlc(t,pa,pb,U):
    origin_t=t
    D=int(round(np.shape(pa)[0]**(1/3)))
    U1=np.kron(U,np.eye(D))
    U2=np.kron(np.eye(D),U)
    pb=U1.dot(pb.dot(U1.conj().transpose()))
    pbnew=np.zeros((D**2,D**2),dtype=complex)
    zero_vector=np.zeros(D,dtype=complex)
    for i in range(D):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(D**2))
        pbnew += project.dot(pb.dot(project.conj().transpose()))
    pb=pbnew
    t=t-1
    while t>0:
        pb=lightconechannelthreesite(U,pb)
        t=t-1
    pb=np.kron(pb,np.eye(D))
    pb=U2.dot(pb.dot(U2.conj().transpose()))
    return np.trace(pa.dot(pb))/D**(origin_t+0.0)

def cf3timechannelleft(U,a):
    dim=int(round(np.shape(a)[0]**(1/3)))
    U=np.kron(U,U)
    a=np.kron(a,np.eye(dim))
    a=U.dot(a.dot(U.conj().T))
    return_matrix=np.zeros((dim**3,dim**3),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim**3),new_vector)
        return_matrix += project.dot(a.dot(project.conj().transpose()))
    return return_matrix

def cf3timechannelright(U,a):
    dim=int(round(np.shape(a)[0]**(1/3)))
    U=np.kron(U,U)
    a=np.kron(np.eye(dim),a)
    a=U.dot(a.dot(U.conj().T))
    return_matrix=np.zeros((dim**3,dim**3),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim**3))
        return_matrix += project.dot(a.dot(project.conj().transpose()))
    return return_matrix

def cf3sitesleft1(t,pa,pb,U):
    origin_t=t
    D=int(round(np.shape(pa)[0]**(1/3)))
    newpb=np.zeros([D**2,D**2],dtype=complex)
    zero_vector=np.zeros(D,dtype=complex)
    for i in range(D):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(D**2),new_vector)
        newpb += project.dot(pb.dot(project.conj().transpose()))
    pb=newpb
    pb=U.dot(pb.dot(U.conj().transpose()))
    t=t-1
    while(t>1):
        pb=channeltwosite(U,pb)
        pb=U.dot(pb.dot(U.conj().transpose()))
        t=t-2
    if t==0:
        ans=np.trace(pa.dot(np.kron(np.eye(D),pb)))
        return ans/(D**(origin_t+0.0))
    if t==1:
        pb=np.kron(np.eye(D),pb)
        pb=np.kron(pb,np.eye(D))
        Unew=np.kron(U,U)
        pb=Unew.dot(pb.dot(Unew.conj().transpose()))
        pa=np.kron(pa,np.eye(D))
        ans=np.trace(pa.dot(pb))
        return ans/(D**(origin_t+0.0))


def cf3sitesright1(t,pa,pb,U):
    pb=cf3timechannelleft(U,pb)
    origin_t=t
    t=t-1
    D=int(round(np.shape(pa)[0]**(1/3)))
    newpb=np.zeros([D**2,D**2],dtype=complex)
    zero_vector=np.zeros(D,dtype=complex)
    for i in range(D):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(D**2))
        newpb += project.dot(pb.dot(project.conj().transpose()))
    pb=newpb
    if t>0:
        pb=U.dot(pb.dot(U.conj().transpose()))
        t=t-1
    while t>1:
        pb=channeltwosite(U,pb)
        pb=U.dot(pb.dot(U.conj().T))
        t=t-2
    if t==0:
        pb=np.kron(pb,np.eye(D))
        ans=np.trace(pa.dot(pb))
        return ans/(D**(origin_t+0.0))
    if t==1:
        pb=np.kron(pb,np.eye(D))
        pb=cf3timechannelright(U,pb)
        ans=np.trace(pa.dot(pb))
        return ans/(D**(origin_t+0.0))



def cf3sitestime(t,pa,pb,U):
    dim=int(round(pa.shape[0]**(1/3)))
    origin_t=t
    while t>0:
        pb=cf3timechannelleft(U,pb)
        t=t-1
        if t>0:
            pb=cf3timechannelright(U,pb)
            t=t-1
    return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def cf3sitesleft2(t,pa,pb,U):
    origin_t=t
    dim=int(round(np.shape(pa)[0]**(1/3)))
    if t%2==0:
        zero_vector=np.zeros(dim,dtype=complex)
        newpb=np.zeros([dim**2,dim**2],dtype=complex)
        for i in range(dim):
            new_vector=np.copy(zero_vector)
            new_vector[i]=1
            project=np.kron(np.eye(dim**2),new_vector)
            newpb += project.dot(pb.dot(project.conj().transpose()))
        pb=newpb
        newpa=np.zeros([dim**2,dim**2],dtype=complex)
        for i in range(dim):
            new_vector=np.copy(zero_vector)
            new_vector[i]=1
            project=np.kron(new_vector,np.eye(dim**2))
            newpa += project.dot(pa.dot(project.conj().transpose()))
        pa=newpa
        return twositeoperator1leftstep(t,pa,pb,U)/dim
    zero_vector=np.zeros(dim,dtype=complex)
    newpb=np.zeros([dim**2,dim**2],dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim**2),new_vector)
        newpb += project.dot(pb.dot(project.conj().transpose()))
    pb=newpb
    pb=U.dot(pb.dot(U.conj().transpose()))
    t=t-1
    newpb=np.zeros([dim,dim],dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(np.eye(dim),new_vector)
        newpb+=project.dot(pb.dot(project.conj().transpose()))
    pb=newpb
    if t==0:
        return np.trace(pa.dot(np.kron(np.eye(dim**2),pb)))/(dim**2)
    while t>2:
        pb=epsilonL(U,pb)
        pb=epsilonR(U,pb)
        t=t-2
    pb=U.dot(np.kron(np.eye(dim),pb).dot(U.conj().transpose()))
    pb=np.kron(np.eye(dim),pb)
    pb=cf3timechannelleft(U,pb)
    return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def cf3sitesright2(t,pa,pb,U):
    dim=int(round(np.shape(pa)[0]**(1/3)))
    origin_t=t
    if t%2==1:
        if t==1:
            pb=np.kron(pb,np.eye(dim))
            pb=np.kron(U,U).dot(pb.dot(np.kron(U,U).conj().transpose()))
            pb=np.kron(pb,np.eye(dim))
            pa=np.kron(np.eye(dim**2),pa)
            return np.trace(pa.dot(pb))/dim**2
        pb=cf3timechannelleft(U,pb)
        t=t-1
        zero_vector=np.zeros(dim,dtype=complex)
        newpb=np.zeros([dim**2,dim**2],dtype=complex)
        for i in range(dim):
            new_vector=np.copy(zero_vector)
            new_vector[i]=1
            project=np.kron(new_vector,np.eye(dim**2))
            newpb += project.dot(pb.dot(project.conj().transpose()))
        pb=newpb
        newpa=np.zeros([dim**2,dim**2],dtype=complex)
        for i in range(dim):
            new_vector=np.copy(zero_vector)
            new_vector[i]=1
            project=np.kron(np.eye(dim**2),new_vector)
            newpa += project.dot(pa.dot(project.conj().transpose()))
        pa=newpa
        return twositeoperator1rightstep(t,pa,pb,U)/dim**2
    pb=cf3timechannelleft(U,pb)
    t=t-1
    zero_vector=np.zeros(dim,dtype=complex)
    newpb=np.zeros([dim**2,dim**2],dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim**2))
        newpb += project.dot(pb.dot(project.conj().transpose()))
    pb=newpb
    pb=U.dot(pb.dot(U.conj().transpose()))
    t=t-1
    newpb=np.zeros([dim,dim],dtype=complex)
    for i in range(dim):
        new_vector=np.copy(zero_vector)
        new_vector[i]=1
        project=np.kron(new_vector,np.eye(dim))
        newpb+=project.dot(pb.dot(project.conj().transpose()))
    pb=newpb
    while t>2:
        pb=epsilonR(U,pb)
        pb=epsilonL(U,pb)
        t=t-2
    pb=np.kron(pb,np.eye(dim)).T
    pb=U.dot(pb.dot(U.conj().T))
    pb=np.kron(pb,np.eye(dim))
    pb=cf3timechannelright(U,pb)
    return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

    



def support3correlationfunctionqubitcase(t,j,i,pa,pb,U):
    if i%2!=0:
        print("i must be even")
        return
    if j==i+t:
        return cf3sitesrightlc(t,pa,pb,U)
    if j==i:
        return cf3sitestime(t,pa,pb,U)
    if j==i+1:
        return cf3sitesright1(t,pa,pb,U)
    if j==i-1:
        return cf3sitesleft1(t,pa,pb,U)
    if j==i-2:
        return cf3sitesleft2(t,pa,pb,U)
    if j==i+2:
        return cf3sitesright2(t,pa,pb,U)
    return 0



def correlation_function_on_lightcone(t,L,pa,pb,U1,U2):
    dim=np.shape(pa)[0]
    origin_t=t
    if L==0:
        while t>0:
            pb=ML(U1,pb)
            t=t-1
            if t>0:
                pb=ML(U2,pb)
                t=t-1
    else:
        while t>0:
            pb=MR(U1,pb)
            t=t-1
            if t>0:
                pb=MR(U2,pb)
                t=t-1
    return np.trace(pa.dot(pb))/(dim**(origin_t+0.0))

def uniform_matrix(t,L,i,j,U1,U2):
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    origin_t=t
    pauli_group=[X,Y,Z]
    pa=pauli_group[i]
    pb=pauli_group[j]
    if L==0:
        while t>0:
            pb=epsilonR(U1,pb)
            t=t-1
            if t>0:
                pb=epsilonL(U2,pb)
                t=t-1
    else:
        while t>0:
            pb=epsilonL(U2,pb)
            t=t-1
            if t>0:
                pb=epsilonR(U1,pb)
                t=t-1
    return np.trace(pa.dot(pb))/(2.0**(origin_t+1.0))


def translation_operator(d, N, k=1):
    """
    Generates the translation operator for a quantum system with local dimension `d` and `N` sites.
    
    Parameters:
    d (int): the local dimension of the Hilbert space.
    N (int): the total number of sites in the system.
    k (int): the number of sites to translate.
    
    Returns:
    np.ndarray: the translation operator.
    it translation i000000 to 000000i
    """
    
    # Create the translation operator
    T = np.zeros((d**N, d**N), dtype=np.complex128)
    for i in range(d**N):
        # Compute the translation index for the i-th basis state
        ti = (i // (d**(N-k))) + ((i % (d**(N-k))) * d**(k))
        T[ti, i] = 1.0
    
    return T


class correlation_function_direct:
    def __init__(self,dim,length):
        self.dim=dim
        self.total_matrix=None
        self.first_layer=None
        self.second_layer=None
        self.uniform_matrix=None
        self.length=length
        self.TO=translation_operator(dim,length)

    def calculate_correlation_function_direct(self,pa,pb,U=None):
        if U is None:
            U=self.total_matrix
        return np.trace(pa.dot(U.dot(pb.dot(U.conj().transpose()))))/(self.dim**self.length)


    def create_layer_first(self,U=None):
        if U is None:
            U=self.uniform_matrix
        self.first_layer=np.eye(1)
        for i in range(self.length//2):
            self.first_layer=np.kron(self.first_layer,U)
    
    def create_second_layer(self,U=None):
        if U is None:
            U=self.uniform_matrix
        self.second_layer=np.eye(1)
        for i in range(self.length//2):
            self.second_layer=np.kron(self.second_layer,U)
        self.second_layer=self.TO.dot(self.second_layer.dot(self.TO.conj().transpose()))

    def create_totol_matrix(self, t, first_layer=None, second_layer=None):
        if first_layer is None:
            first_layer=self.first_layer
        if second_layer is None:
            second_layer=self.second_layer
        self.total_matrix=np.eye(self.dim**self.length)
        time=t
        while time>0:
            self.total_matrix=first_layer.dot(self.total_matrix)
            time=time-1
            if time>0:
                self.total_matrix=second_layer.dot(self.total_matrix)
                time=time-1
    
    def create_uniform_matrix(self,U):
        self.uniform_matrix=U
    

def self_kron(oplist):
    returnop=oplist[0]
    for i in range(1,len(oplist)):
        returnop=np.kron(returnop,oplist[i])
    return returnop

def general_extend_matrix(M,D,L,i):
    translationmatrix=translation_operator(D,L,i)
    Md=np.shape(M)[0]
    returnmatrix=np.kron(M,np.eye(D**L//Md))
    return translationmatrix.conj().transpose().dot(returnmatrix.dot(translationmatrix))

def extend_matrix(M, D, L, i):
    """
    Extends a `D x D` matrix `M` to a `D**L x D**L` matrix that acts on the `i`-th site of a `D`-dimensional
    Hilbert space with `L` sites.

    Parameters:
    M (np.ndarray): the matrix to extend.
    D (int): the local dimension of the Hilbert space.
    L (int): the total number of sites in the system.
    i (int): the index of the site to act on (0 <= i < L).
    
    Returns:
    np.ndarray: the extended matrix.
    """
    # Create a list of D x D identity matrices for each site
    Id_list = [np.eye(D) for _ in range(L)]
    
    # Replace the identity matrix at the i-th site with M
    Id_list[i] = M
    
    # Compute the tensor product of the identity matrices
    T = self_kron(Id_list)
    
    return T


def correlation_function_uniform_2ndH(U1,U2,t,pa,pb,lengthL,initialsite,fp=None):
    pb=extend_matrix(pb,2,lengthL,initialsite)
    intm=interaction(0,0,np.pi/4)
    U=np.kron(U1,U2).dot(intm)
    CFC=correlation_function_direct(2,lengthL)
    CFC.create_uniform_matrix(U)
    CFC.create_layer_first()
    CFC.create_second_layer()
    CFC.create_totol_matrix(t)
    for i in range(lengthL):
        panew=extend_matrix(pa,2,lengthL,i)
        resultans=CFC.calculate_correlation_function_direct(panew,pb)
        print("{}\t{}\n".format(i,resultans))
        if fp is not None:
            fp.write("{}\t{}\t{}\n".format(i,t,resultans))

def correlation_function_square_1site(U,t,lengthL,initialsite,fp=None):
    h=np.random.rand(2,2)+1j*np.random.rand(2,2)
    h=(h+h.conj().T)/2
    tr=np.trace(h)
    h=h-tr*np.eye(2)/2
    pb=h
    pa=h
    pb=extend_matrix(pb,2,lengthL,initialsite)
    CFC=correlation_function_direct(2,lengthL)
    CFC.create_uniform_matrix(U)
    CFC.create_layer_first()
    CFC.create_second_layer()
    CFC.create_totol_matrix(t)
    for i in range(lengthL):
        panew=extend_matrix(pa,2,lengthL,i)
        resultans=CFC.calculate_correlation_function_direct(panew,pb)
        print("{}\t{}\n".format(i,resultans))
        if fp is not None:
            fp.write("{}\t{}\t{}\n".format(i,t,resultans))

def frequentlyused2nd():
    intm=interaction(0,0,np.pi/4)
    r=1.24056
    theta=np.arcsin(1/(np.sqrt(2)*np.sin(r)))
    phi=-0.4764
    U1=singleoperator(r,theta,phi)
    U2=singleoperator(r,theta,phi)
    U=np.kron(U1,U2).dot(intm)
    return U

def another2nd():
    intm=interaction(0,0,np.pi/4)
    r=np.pi/4
    theta=np.pi/2
    phi=0
    U1=singleoperator(r,theta,phi)
    U2=singleoperator(r,theta,phi)
    U=np.kron(U1,U2).dot(intm)
    return U

def Random2_2nd():
    intm=interaction(0,0,np.pi/4)
    r=1
    theta=np.random.rand()*np.pi
    phi=3
    U1=singleoperator(r,theta,phi)
    U2=singleoperator(r,theta,phi)
    U=np.kron(U1,U2).dot(intm)
    return U


def correlation_function_2sites(U1,U2,t,pa,pb,lengthL,initialsite,fp=None):
    pb=extend_matrix(pb,4,int(lengthL/2),initialsite)
    intm=interaction(0,0,np.pi/4)
    U=np.kron(U1,U2).dot(intm)
    CFC=correlation_function_direct(2,lengthL)
    CFC.create_uniform_matrix(U)
    CFC.create_layer_first()
    CFC.create_second_layer()
    CFC.create_totol_matrix(t)
    for i in range(int(lengthL/2)):
        panew=extend_matrix(pa,4,int(lengthL/2),i)
        resultans=CFC.calculate_correlation_function_direct(panew,pb)
        print("{}\t{}\n".format(i,resultans))
        if fp is not None:
            fp.write("{}\t{}\t{}\n".format(i,t,resultans))






def test_zhiyuanwangidea(t,L,h=None):
    if h is None:
        h=np.random.rand(2,2)+1j*np.random.rand(2,2)
        h=(h+h.conj().T)/2
        tr=np.trace(h)
        h=h-tr*np.eye(2)/2
    pa=h
    pb=h
    pb=extend_matrix(pb,2,L,4)
    U=interaction(0,-3*np.pi/4,-np.pi/2)
    u1=singleoperator(-1.47219,-1.26402,-3.7839)
    u2=singleoperator(2.15135,-2.17328,5.57852)
    U=np.kron(u1,u2).dot(U)
    CFC=correlation_function_direct(2,L)
    CFC.create_uniform_matrix(U)
    CFC.create_layer_first()
    CFC.create_second_layer()
    CFC.create_totol_matrix(t)
    for i in range(L):
        panew=extend_matrix(pa,2,L,i)
        resultans=CFC.calculate_correlation_function_direct(panew,pb)
        print("{}\t{}\n".format(i,resultans))




def correlation_function_multisites(U1,U2,t,pa,pb,lengthL,initialsite,fp=None):
    pb=general_extend_matrix(pb,2,lengthL,initialsite)
    intm=interaction(0,0,np.pi/4)
    U=np.kron(U1,U2).dot(intm)
    CFC=correlation_function_direct(2,lengthL)
    CFC.create_uniform_matrix(U)
    CFC.create_layer_first()
    CFC.create_second_layer()
    CFC.create_totol_matrix(t)
    for i in range(lengthL):
        panew=general_extend_matrix(pa,2,lengthL,i)
        resultans=CFC.calculate_correlation_function_direct(panew,pb)
        print("{}\t{}\n".format(i,resultans))
        if fp is not None:
            fp.write("{}\t{}\t{}\n".format(i,t,resultans))

def random_matrixU2(t,pa,pb,length,ini):
    r=1.24056
    theta=np.arcsin(1/(np.sqrt(2)*np.sin(r)))
    phi=-0.4764
    U1=singleoperator(r,theta,phi)
    U2=unitary_group.rvs(2)
    correlation_function_uniform_2ndH(U1,U2,t,pa,pb,length,ini)

def easyway_t_2sites(t,L,sites):
    h = np.random.rand(4,4) + 1j*np.random.rand(4,4)
    h = (h + h.conj().T) / 2
    tr = np.trace(h)
    h = h - tr*np.eye(4)/4
    r=1.24056
    theta=np.arcsin(1/(np.sqrt(2)*np.sin(r)))
    phi=-0.4764
    u=singleoperator(r,theta,phi)
    correlation_function_2sites(u,u,t,h,h,L,sites)

def easyway_t_3sites(t,L,initialsite):
    h=np.random.rand(8,8)+1j*np.random.rand(8,8)
    h=(h+h.conj().T)/2
    tr=np.trace(h)
    h=h-tr*np.eye(8)/8
    r=1.24056
    theta=np.arcsin(1/(np.sqrt(2)*np.sin(r)))
    phi=-0.4764
    u=singleoperator(r,theta,phi)
    correlation_function_multisites(u,u,t,h,h,L,initialsite)

def decay_of_correlationfunction():
    h=np.load("./plot_data/correlation_function_support3.npy")
    fp=open("./plot_data/decay_of_correlationfunction.txt",'w')
    U=frequentlyused2nd()
    for t in range(1,100):
        ans1=support3correlationfunctionqubitcase(t,20,20,h,h,U)
        ans2=support3correlationfunctionqubitcase(t,20+t,20,h,h,U)
        fp.write("{}\t{}\t{}\n".format(t,abs(ans1),abs(ans2)))
    fp.close()


def create_picture_datasupport3(filename):
    h=np.random.rand(8,8)+1j*np.random.rand(8,8)
    h=(h+h.conj().T)/2
    tr=np.trace(h)
    h=h-tr*np.eye(8)/8
    normalization=np.trace(h.conj().T.dot(h))
    h=h/np.sqrt(normalization)
    U=frequentlyused2nd()
    np.save("./plot_data/"+filename.replace(".txt",".npy"),h)
    fp=open("./plot_data/"+filename,'w')
    pa=h
    pb=h
    for t in range(1,21):
        for j in range(0,41):
            tempans=support3correlationfunctionqubitcase(t,j,20,pa,pb,U)
            fp.write("{}\t{}\t{}\n".format(t,j,tempans))
            print("{}\t{}\t{}\n".format(t,j,tempans))
    fp.close()

def create_picture_datasupport2(filename):
    U=generate_matrix(6,0,3,0)
    U1=np.zeros([3,3],dtype=complex)
    X=np.array([[0,1],[1,0]])
    for i in range(3):
        U1[i,i]=np.exp(1j*2*np.pi*np.random.random())
    U2=np.zeros([3,3],dtype=complex)
    for i in range(3):
        U2[i,i]=np.exp(1j*2*np.pi*np.random.random())
    U1=np.kron(np.eye(2),U1)
    U2=np.kron(np.eye(2),U2)
    U1=np.kron(X,unitary_group.rvs(3))
    U2=np.kron(X,unitary_group.rvs(3))
    U1=uniform_group.rvs(6)
    U2=uniform_group.rvs(6)
    U=np.kron(U1,U2).dot(U)
    a,b=generalwhether_left_invariant(U)
    print("the left deviation is {}".format(np.trace(b.conj().T.dot(b))))
    a,b=generalwhether_right_invariant(U)
    print("the right deviation is {}".format(np.trace(b.conj().T.dot(b))))
    h=np.random.rand(36,36)+1j*np.random.rand(36,36)
    h=(h+h.conj().T)/2
    tr=np.trace(h)
    h=h-tr*np.eye(36)/36
    normalization=np.trace(h.conj().T.dot(h))
    h=h/np.sqrt(normalization)
    np.save("./plot_data/"+filename.replace(".txt",".npy"),h)
    fp=open("./plot_data/"+filename,'w')
    pa=h
    pb=h
    ts=[]
    xs=[[],[],[]]
    tsfit=[]
    xsfit=[]
    for t in range(20,60,2):
        tsfit.append(t)
        tempans=support2correlationfunction(t,20,20,pa,pb,U)
        xsfit.append(np.abs(tempans))
    plt.plot(tsfit,xsfit,label="fit")
    plt.show()
    tsfit=np.log(tsfit)
    xsfit=np.log(xsfit)
    plt.plot(tsfit,xsfit,label="fit")
    plt.show()
    for t in range(1,21):
        for j in range(0,41):
            tempans=support2correlationfunction(t,j,20,pa,pb,U)
            fp.write("{}\t{}\t{}\n".format(t,j,tempans))
            if j==19:
                print("{}\t{}\t{}\n".format(t,j,tempans))
                ts.append(t)
                xs[0].append(np.abs(tempans))
            if j==20:
                print("{}\t{}\t{}\n".format(t,j,tempans))
                xs[1].append(np.abs(tempans))
            if j==21:
                print("{}\t{}\t{}\n".format(t,j,tempans))
                xs[2].append(np.abs(tempans))
    plt.plot(ts,xs[0],label="x=19")
    plt.plot(ts,xs[1],label="x=20")
    plt.plot(ts,xs[2],label="x=21")
    plt.legend()
    plt.show()
    fp.close()


def create_picture_data():
    r=1.24056
    theta=np.arcsin(1/(np.sqrt(2)*np.sin(r)))
    phi=-0.4764
    U1=singleoperator(r,theta,phi)
    U2=unitary_group.rvs(2)
    print(U2)
    np.savetxt("./plot_data/correlation_function_randomU2.txt",U2,delimiter="\t")
    fp=open("./plot_data/correlation_function_randomU2.txt","a")
    X=np.array([[0,1],[1,0]])
    Z=np.array([[1,0],[0,-1]])
    pb=X
    pa=1000*X+1000*Z
    for t in range(0,6):
        correlation_function_uniform_2ndH(U1,U2,t,pa,pb,14,5,fp)
    fp.close()

def compare_channel_easy3(t,L,initial):
    h=np.random.rand(8,8)+1j*np.random.rand(8,8)
    h=(h+h.conj().T)/2
    tr=np.trace(h)
    h=h-tr*np.eye(8)/8
    r=1.24056
    theta=np.arcsin(1/(np.sqrt(2)*np.sin(r)))
    phi=-0.4764
    u=singleoperator(r,theta,phi)
    correlation_function_multisites(u,u,t,h,h,L,initial)
    U=frequentlyused2nd()
    print("the channel result is")
    for j in range(L):
        tempans=support3correlationfunctionqubitcase(t,j,initial,h,h,U)/8.0
        print("{}\t{}\n".format(j,tempans))

def quantum_quench_from_initial_state2site(statematrix,h,U,t):
    if t==0:
        return np.trace(statematrix.dot(h))
    statematrix=channeltwosite(U,statematrix)
    t=t-1
    return twositeoperator(t,h,statematrix,U,U)

def thermalize_quantum_quench(filename,h=None,U=None):
    statematrix=np.array([[0,0,0,0],[0,1,1,0],[0,1,1,0],[0,0,0,0]])/2
    if h is None:
        h=np.random.rand(4,4)+1j*np.random.rand(4,4)
        h=(h+h.conj().T)/2
        tr=np.trace(h)
        h=h-tr*np.eye(4)/4
        norm=np.trace(h.conj().T.dot(h))
        h=h/np.sqrt(norm)
    if U is None:
        U=frequentlyused2nd()
    np.save("./plot_data/"+filename.replace(".txt",".npy"),h)
    fp=open("./plot_data/"+filename,'w')
    for t in range(0,61):
        ans=quantum_quench_from_initial_state2site(statematrix,h,U,t)
        fp.write("{}\t{}\n".format(t,ans))
        print("{}\t{}\n".format(t,ans))
    fp.close()



if __name__=='__main__':
    if sys.argv[1]=='1':
        create_picture_datasupport3("correlation_function_support3.txt")
        decay_of_correlationfunction
    if sys.argv[1]=='2':
        create_picture_datasupport2("correlation_function_support2trialfile.txt")
    if sys.argv[1]=='3':
        h=np.random.rand(4,4)+1j*np.random.rand(4,4)
        h=(h+h.conj().T)/2
        tr=np.trace(h)
        h=h-tr*np.eye(4)/4
        norm=np.trace(h.conj().T.dot(h))
        h=h/np.sqrt(norm)
        print("Random1")
        U=frequentlyused2nd()
        thermalize_quantum_quench("thermalize_quantum_quench.txt",h,U)
        print("Non Ergodic")
        U=another2nd()
        thermalize_quantum_quench("thermalize_quantum_quenchanother.txt",h,U)
        print("Random2")
        U=Random2_2nd()
        if h is not None and U is not None:
            thermalize_quantum_quench("thermalize_quantum_quenchrandom2.txt",h,U)

        else:
            print("Error: h or U is None")
