from qutip import *
import numpy as np
import matplotlib.pyplot as plt


def PD(i,j,N):
    return min(abs(i-j),N-abs(i-j))


def LongrangeHamiltonian(N,alpha):
    normfactor=np.max([N**(1-alpha),1])
    X = sigmax()
    Y = sigmay()
    Z = sigmaz()
    returnop=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        for j in range(i+1,N):
            returnop+=1.0*tensor([identity(2)]*i+[X]+[identity(2)]*(j-i-1)+[X]+[identity(2)]*(N-j-1))/normfactor/((1.0+PD(i,j,N))**alpha)
            returnop+=1.0*tensor([identity(2)]*i+[Y]+[identity(2)]*(j-i-1)+[Y]+[identity(2)]*(N-j-1))/normfactor/((1.0+PD(i,j,N))**alpha)
            returnop+=1.0*tensor([identity(2)]*i+[Z]+[identity(2)]*(j-i-1)+[Z]+[identity(2)]*(N-j-1))/normfactor/((1.0+PD(i,j,N))**alpha)
    return returnop

def NearestneighbourHamiltonian(N):
    X = sigmax()
    returnop=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N-1):
        returnop+=1.0*tensor([identity(2)]*i+[X]+[X]+[identity(2)]*(N-i-2))
    returnop+=-1.0*tensor([X]+[identity(2)]*(N-2)+[X])
    return returnop

def groundstate(H):
    states=H.groundstate()
    print("The ground state energy is {}".format(states[0]))
    hardmard_transform=Qobj(dims=[[2],[2]])
    hardmard_transform.data[0,0]=1.0/np.sqrt(2)
    hardmard_transform.data[0,1]=1.0/np.sqrt(2)
    hardmard_transform.data[1,0]=1.0/np.sqrt(2)
    hardmard_transform.data[1,1]=-1.0/np.sqrt(2)
    hardmard=tensor([hardmard_transform]*10)
    print("The ground state is {}".format(hardmard*states[1]))
    t1=((hardmard*states[1]).data[int('1010101010',2),0])
    t2=((hardmard*states[1]).data[int('0101010101',2),0])
    print(abs(t1)**2+abs(t2)**2)
    return states[1]

def decay_correlationfunction(psi,N):
    X=sigmax()
    distances=[]
    correlators=[]
    for i in range(1,N):
        distances.append(i)
        correlators.append(expect(tensor([X]+[identity(2)]*(i-1)+[X]+[identity(2)]*(N-i-1)),psi)*(-1)**i)
    plt.plot(distances,correlators)
    plt.show()
    print("The single site expectation value is {},{}".format(expect(tensor([X]+[identity(2)]*(N-1)),psi),expect(tensor([identity(2)]*(N-1)+[X]),psi)))


def different_alpha_correlators(N):
    for alpha in range(21):
        print(alpha/10.0)
        sum_energy=0
        for i in range(1,N):
            sum_energy+=0.5*((-1)**(i))*(N**(alpha/10.0))/(PD(i,0,N)**(alpha/10.0))
        print("The naive ground state energy is {}".format(sum_energy))
        psi=groundstate(LongrangeHamiltonian(N,alpha/10.0))
        decay_correlationfunction(psi,N)

if __name__=='__main__':
    different_alpha_correlators(10)

