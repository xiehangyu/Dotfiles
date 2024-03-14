import numpy as np
import numpy.linalg as LA
from ncon import ncon
from scipy.stats import unitary_group

def one_iteration_of_U(tensorABC,U_origin):
    #The input tensor should be 3-ranks with ABC. The middle B is the one to be splitted
    #For the returnned U, the input is two legs and the output is one leg
    dA, dB, dC = tensorABC.shape
    mediate_matrix=ncon([tensorABC,tensorABC.conj(),U_origin.conj(),U_origin,tensorABC,tensorABC.conj(),U_origin.conj()],[[1,-1,2],[3,4,2],[4,5,-3],[7,5,6],[3,7,8],[1,9,8],[9,-2,6]])
    u,s,v=LA.svd(mediate_matrix.reshape(dB,dB))
    U_new=u.conj()@v.conj()
    return U_new.reshape([dB,int(np.sqrt(dB)),int(np.sqrt(dB))])

def purity_calculation(tensorABC,U_origin):
    return ncon([tensorABC,tensorABC.conj(),U_origin.conj(),U_origin,tensorABC,tensorABC.conj(),U_origin.conj(),U_origin],[[1,10,2],[3,4,2],[4,5,12],[7,5,6],[3,7,8],[1,9,8],[9,11,6],[10,11,12]])

def BestU(tensorABC):
    epsilon=1E-6
    norm=LA.norm(tensorABC)
    tensorABC=tensorABC/norm
    dA,dB,dC=tensorABC.shape
    U=unitary_group.rvs(dB).reshape([dB,int(np.sqrt(dB)),int(np.sqrt(dB))])
    #U=np.eye(dB).reshape([dB,int(np.sqrt(dB)),int(np.sqrt(dB))])
    purity1=purity_calculation(tensorABC,U)
    print("The original purity is{}".format(purity1))
    purity2=-100
    while abs(purity1-purity2)>epsilon:
        U=one_iteration_of_U(tensorABC,U)
        purity2=purity1
        purity1=purity_calculation(tensorABC,U)
        print("{}".format(purity1))
    print("The purity is{}".format(purity1))
    tensorABC=tensorABC*norm
    return U

def addtriangle(tensorABC,chi):
    '''
    B           B
        C to       triangle   C

    A           A
    '''
    dA,dB,dC=tensorABC.shape
    u,s,v=LA.svd(tensorABC.reshape([dA,dB*dC]))
    print(s)
    tensorA=u[:,:chi**2]
    tensorBC=(np.diag(s[:chi**2])@v[:chi**2,:]).reshape([chi**2,dB,dC])
    temp_tensor=tensorBC.transpose([1,0,2])
    U=BestU(temp_tensor)
    tensorA=ncon([tensorA,U.conj()],[[-1,1],[1,-2,-3]])
    tensorBC=ncon([U,tensorBC],[[1,-1,-2],[1,-3,-4]])
    u,s,v=LA.svd(tensorBC.transpose([0,2,1,3]).reshape([chi*dB,chi*dC]))
    print(s)
    tensorB=(u[:,:chi]@np.diag(s[:chi])).reshape([chi,dB,chi]).transpose([1,0,2])
    tensorC=v[:chi,:].reshape([chi,chi,dC]).transpose([2,0,1])
    return tensorA, tensorB, tensorC

def addtrianglewithoutboundary(tensorABC,chi):
    '''
    B           B
        C to       triangle   C

    A           A
    '''
    dA,dB,dC=tensorABC.shape
    u,s,v=LA.svd(tensorABC.reshape([dA,dB*dC]))
    print(s)
    cutoff1=min(int(np.sqrt(dB*dC)),int(np.sqrt(dA)))
    tensorA=u[:,:cutoff1**2]
    tensorBC=(np.diag(s[:cutoff1**2])@v[:cutoff1**2,:]).reshape([cutoff1**2,dB,dC])
    temp_tensor=tensorBC.transpose([1,0,2])
    U=BestU(temp_tensor)
    tensorA=ncon([tensorA,U.conj()],[[-1,1],[1,-2,-3]])
    tensorBC=ncon([U,tensorBC],[[1,-1,-2],[1,-3,-4]])
    u,s,v=LA.svd(tensorBC.transpose([0,2,1,3]).reshape([cutoff1*dB,cutoff1*dC]))
    print(s)
    cutoff2=min([cutoff1*dB,cutoff1*dC])
    tensorB=(u[:,:cutoff2]@np.diag(s[:cutoff2])).reshape([cutoff1,dB,cutoff2]).transpose([1,0,2])
    tensorC=v[:chi,:].reshape([chi,chi,dC]).transpose([2,0,1])
    return tensorA, tensorB, tensorC


def Moses_Move(C_arrays,chi,D):
    Lx=len(C_arrays)
    bottom_lay_C_arrays=[]
    Upper_lay_C_arrays=[]
    tensorA,tensorB,tensorC=addtriangle(C_arrays[0].transpose([0,1,2,4,3]).reshape(D*chi**2,chi,chi),chi)
    bottom_lay_C_arrays.append(tensorA.reshape([D,chi,chi,chi,chi]).transpose([0,1,2,4,3]))
    Upper_lay_C_arrays.append(tensorB.transpose([1,2,0]))
    newtensor=(ncon([tensorC,C_arrays[1]],[[1,-4,-2],[-1,1,-3,-6,-5]])).reshape([D*chi**2,chi**2,chi])
    for i in range(1,Lx-1):
        tensorA,tensorB,tensorC=addtriangle(newtensor,chi)
        bottom_lay_C_arrays.append(tensorA.reshape([D,chi,chi,chi,chi]).transpose([0,1,2,4,3]))
        Upper_lay_C_arrays.append(tensorB.reshape([chi,chi,chi,chi]).transpose([0,2,3,1]))
        newtensor=(ncon([tensorC,C_arrays[i+1]],[[1,-4,-2],[-1,1,-3,-6,-5]])).reshape([D*chi**2,chi**2,chi])
    u,s,v=LA.svd(newtensor.transpose([0,2,1]).reshape([D*chi**3,chi**2]))
    bottom_lay_C_arrays.append((u[:,:chi]).reshape([D,chi,chi,chi,chi]))
    Upper_lay_C_arrays.append((np.diag(s[:chi])@v[:chi,:]).reshape([chi,chi,chi]).transpose([1,0,2]))
    return bottom_lay_C_arrays,Upper_lay_C_arrays

    
def test_Moses_Move(C_arrays,chi,D):
    original_tensor=C_arrays[0]
    Lx=len(C_arrays)
    for i in range(1,Lx):
        original_tensor=ncon([original_tensor,C_arrays[i]],[[-1]+[-2-j for j in range(2*i)]+[1]+[-2*i-6-j for j in range(i)],[-2*i-2,1,-2*i-3,-2*i-4,-2*i-5]])
    bottom_lay_C_arrays,Upper_lay_C_arrays=Moses_Move(C_arrays,chi,D)
    new_tensor=ncon([bottom_lay_C_arrays[0],Upper_lay_C_arrays[0]],[[-1,-2,-3,-4,1],[1,-5,-6]])
    for i in range(1,Lx-1):
        new_tensor=ncon([new_tensor,bottom_lay_C_arrays[i],Upper_lay_C_arrays[i]],[[-1]+[-2-j for j in range(2*i)]+[1,2]+[-2*i-7-j for j in range(i)],[-2*i-2,1,-2*i-3,-2*i-4,3],[2,3,-2*i-5,-2*i-6]])
    new_tensor=ncon([new_tensor,bottom_lay_C_arrays[Lx-1],Upper_lay_C_arrays[Lx-1]],[[-1]+[-2-j for j in range(2*Lx-2)]+[1,2]+[-2*Lx-4-j for j in range(Lx-1)],[-2*Lx,1,-2*Lx-1,-2*Lx-2,3],[2,3,-2*Lx-3]])
    origin_norm=LA.norm(original_tensor)
    new_norm=LA.norm(new_tensor)
    difference_norm=LA.norm(original_tensor-new_tensor/new_norm*origin_norm)
    print("The original norm is {}".format(origin_norm))
    print("The new norm is {}".format(new_norm))
    print("The difference norm is {}".format(abs(difference_norm)))
    print("The difference norm ratio is {}".format(abs(difference_norm)/origin_norm/Lx))
    return bottom_lay_C_arrays,Upper_lay_C_arrays

def truncation_of_two_layers(bottom_lay_C_arrays,Upper_lay_C_arrays,chi,D,total_tensor=None):
    if bottom_lay_C_arrays is not None and Upper_lay_C_arrays is not None:
        Lx=len(Upper_lay_C_arrays)
        total_tensor=ncon([Upper_lay_C_arrays[0],bottom_lay_C_arrays[0]],[[-1,-2,1,-5,-6],[-3,-4,1]])
        for i in range(1,Lx-1):
            total_tensor=ncon([total_tensor,Upper_lay_C_arrays[i],bottom_lay_C_arrays[i]],[[-1]+[-2-j for j in range(2*i)]+[2,1]+[-2*i-7-j for j in range(i)],[-2*i-2,1,3,-2*i-5,-2*i-6],[2,-2*i-3,-2*i-4,3]])
        total_tensor=ncon([total_tensor,Upper_lay_C_arrays[Lx-1],bottom_lay_C_arrays[Lx-1]],[[-1]+[-2-j for j in range(2*Lx-2)]+[2,1]+[-2*Lx-4-j for j in range(Lx-1)],[-2*Lx,1,3,-2*Lx-2,-2*Lx-3],[2,-2*Lx-1,3]])
    elif total_tensor is not None:
        Lx=(total_tensor.ndim-2)//3
    else:
        exit(1)
    new_C_arrays=[]
    Lx_new=Lx
    cutoff=chi
    if Lx<=5:
        for i in range(Lx-1):
            firstdim=total_tensor.shape[1]
            enddim=total_tensor.shape[2*Lx_new+1]
            total_tensor=total_tensor.transpose([0,1,2,3*Lx_new+1]+[j for j in range(3,3*Lx_new+1)])
            u,s,v=LA.svd(total_tensor.reshape(D*chi**2*firstdim,D**(Lx_new-1)*chi**(2*Lx_new-2)*enddim))
            u=u[:,:cutoff]
            s=np.diag(np.sqrt(s[:cutoff]))
            v=v[:cutoff,:]
            temp_C_array=(u@s).reshape([D,firstdim,chi,chi,cutoff]).transpose([0,1,2,4,3])
            new_C_arrays.append(temp_C_array)
            Lx_new=Lx_new-1
            total_tensor=(s@v).reshape([cutoff]+[D,chi]*Lx_new+[enddim]+[chi]*(Lx_new)).transpose([1,0]+[i for i in range(2,3*Lx_new+2)])
        new_C_arrays.append(total_tensor)
        return new_C_arrays
    elif Lx%2==1:
        print("The large Lx has to be even")
        exit(1)
    else:
        Lx_half=Lx//2
        firstdim=total_tensor.shape[1]
        enddim=total_tensor.shape[2*Lx+1]
        total_tensor=total_tensor.transpose([i for i in range(2*Lx_half+1)]+[3*Lx+i-Lx_half+2 for i in range(Lx_half)]+[2*Lx_half+1+i for i in range(3*Lx_half+1)])
        u,s,v=LA.svd(total_tensor.reshape([firstdim*(D*chi**2)**Lx_half,(D*chi**2)**Lx_half*enddim]))
        u=u[:,:cutoff]
        s=np.diag(np.sqrt(s[:cutoff]))
        v=v[:cutoff,:]
        total_tensor1=(u@s).reshape([D,firstdim,chi]+[D,chi]*(Lx_half-1)+[chi]*Lx_half+[cutoff]).transpose([i for i in range(2*Lx_half+1)]+[3*Lx_half+1]+[2*Lx_half+1+i for i in range(Lx_half)])
        total_tensor2=(s@v).reshape([cutoff]+[D,chi]*Lx_half+[enddim]+[chi]*Lx_half).transpose([1,0]+[i for i in range(2,3*Lx_half+2)])
        return truncation_of_two_layers(None,None,chi,D,total_tensor1)+truncation_of_two_layers(None,None,chi,D,total_tensor2)


        

def test_truncation_of_two_layers(bottom_lay_C_arrays,Upper_lay_C_arrays,chi,D):
    Lx=len(Upper_lay_C_arrays)
    new_C_arrays=truncation_of_two_layers(bottom_lay_C_arrays,Upper_lay_C_arrays,chi,D)
    new_tensor=new_C_arrays[0]
    for i in range(1,Lx):
        new_tensor=ncon([new_tensor,new_C_arrays[i]],[[-1]+[-2-j for j in range(2*i)]+[1]+[-2*i-6-j for j in range(i)],[-2*i-2,1,-2*i-3,-2*i-4,-2*i-5]])
    new_norm=LA.norm(new_tensor)
    upper_lay_tensor=Upper_lay_C_arrays[0]
    for i in range(1,Lx):
        upper_lay_tensor=ncon([upper_lay_tensor,Upper_lay_C_arrays[i]],[[-1]+[-2-j for j in range(2*i)]+[1]+[-2*i-6-j for j in range(i)],[-2*i-2,1,-2*i-3,-2*i-4,-2*i-5]])
    bottom_lay_tensor=bottom_lay_C_arrays[0]
    for i in range(1,Lx-1):
        bottom_lay_tensor=ncon([bottom_lay_tensor,bottom_lay_C_arrays[i]],[[-j for j in range(1,i+1)]+[1]+[-i-4-j for j in range(i)],[1,-i-1,-i-2,-i-3]])
    bottom_lay_tensor=ncon([bottom_lay_tensor,bottom_lay_C_arrays[Lx-1]],[[-j for j in range(1,Lx)]+[1]+[-Lx-2-j for j in range(Lx-1)],[1,-Lx,-Lx-1]])
    origin_tensor=ncon([bottom_lay_tensor,upper_lay_tensor],[[-3-2*i for i in range(Lx)]+[Lx-i for i in range(Lx)],[-1]+[-i-1 if i%2==1 else i//2 for i in range(1,2*Lx+1)]+[-2*Lx-2-i for i in range(Lx+1)]])
    origin_norm=LA.norm(origin_tensor)
    difference_norm=LA.norm(origin_tensor-new_tensor)
    print("The original norm is {}".format(origin_norm))
    print("The new norm is {}".format(new_norm))
    print("The difference norm is {}".format(abs(difference_norm)))
    print("The difference norm ratio is {}".format(abs(difference_norm)/origin_norm/Lx))
        



if __name__=='__main__':
    Lx=6
    chi=2
    D=2
    C_arrays=[]
    for i in range(Lx):
        C_arrays.append(np.random.random([D,chi,chi,chi,chi])+1j*np.random.random([D,chi,chi,chi,chi]))
    bottom_lay,upper_lay=test_Moses_Move(C_arrays,chi,D)
    test_truncation_of_two_layers(upper_lay,bottom_lay,chi,D)


