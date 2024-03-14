from qutip import *
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
import sys


def PD(i,j,N):
    return min(abs(i-j),N-abs(i-j))
#   return abs(i-j)


def LongrangeHamiltonian(N,alpha):
    normfactor=np.max([N**(1-alpha),1])
    X = sigmax()/2.0
    Y = sigmay()/2.0
    Z = sigmaz()/2.0
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


def site2optimization1D(theta,L,alpha):
    sumans=0
    i=0
    for k in range(1,L):
        if k%2==0:
            sumans += 1/(PD(i,k,L)**alpha)
        else:
            sumans += np.cos(theta)/(PD(i,k,L)**alpha)
    i=1
    for k in range(0,L):
        if k==1:
            continue
        if k%2==0:
            sumans += np.cos(theta)/(PD(i,k,L)**alpha)
        else:
            sumans += 1/(PD(i,k,L)**alpha)
    return sumans

def minimizationsite2D1():
    L=int(sys.argv[1])
    initguess=np.random.random()*2*np.pi
    print("The initguess is {}".format(initguess))
    alpha=float(sys.argv[2])
    print(L,alpha)
    def site2optimization1Dwrapper(theta):
        return site2optimization1D(theta,L,alpha)
    res=minimize(site2optimization1Dwrapper,initguess,method='BFGS',)
    print("最小化的角度是{}".format(res.x))
    print("最小化的值是{}".format(res.fun))
    print(res)
    
def innerproduct_2vector_sphere(theta1,phi1,theta2,phi2):
    return np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2)

def state3optimization1D(thetas,L,alpha):
    theta1=0
    phi1=0
    theta2,phi2,theta3,phi3=thetas
    sumans=0
    i=0
    for k in range(1,L):
        if k%3==0:
            sumans += 1/(PD(i,k,L)**alpha)
        elif k%3==1:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta2,phi2)/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta3,phi3)/(PD(i,k,L)**alpha)
    i=1
    for k in range(L):
        if k==i:
            continue
        if k%3==0:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta2,phi2)/(PD(i,k,L)**alpha)
        elif k%3==1:
            sumans += 1/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta3,phi3)/(PD(i,k,L)**alpha)
    i=2
    for k in range(L):
        if k==i:
            continue
        if k%3==0:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta3,phi3)/(PD(i,k,L)**alpha)
        elif k%3==1:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta3,phi3)/(PD(i,k,L)**alpha)
        else:
            sumans += 1/(PD(i,k,L)**alpha)
    return sumans/3.0*2.0


def minimizationstate3optimization1D():
    L=int(sys.argv[1])
    initguess=np.random.random([4])*2*np.pi
    print("The initguess is {}".format(initguess))
    alpha=float(sys.argv[2])
    print(L,alpha)
    def state3optimization1Dwrapper(thetas):
        return state3optimization1D(thetas,L,alpha)
    res=minimize(state3optimization1Dwrapper,initguess,method='BFGS',)
    print("最小化的角度是{}".format(res.x))
    print("最小化的值是{}".format(res.fun))
    print(res)

    
def state4optimization1D(thetas,L,alpha):
    theta1=0
    phi1=0
    theta2,phi2,theta3,phi3,theta4,phi4=thetas
    sumans=0
    i=0
    for k in range(1,L):
        if k%4==0:
            sumans += 1/(PD(i,k,L)**alpha)
        elif k%4==1:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta2,phi2)/(PD(i,k,L)**alpha)
        elif k%4==2:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta3,phi3)/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta4,phi4)/(PD(i,k,L)**alpha)
    i=1
    for k in range(L):
        if k==i:
            continue
        if k%4==0:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta2,phi2)/(PD(i,k,L)**alpha)
        elif k%4==1:
            sumans += 1/(PD(i,k,L)**alpha)
        elif k%4==2:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta3,phi3)/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta4,phi4)/(PD(i,k,L)**alpha)
    i=2
    for k in range(L):
        if k==i:
            continue
        if k%4==0:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta3,phi3)/(PD(i,k,L)**alpha)
        elif k%4==1:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta3,phi3)/(PD(i,k,L)**alpha)
        elif k%4==2:
            sumans += 1/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta3,phi3,theta4,phi4)/(PD(i,k,L)**alpha)
    i=3
    for k in range(L):
        if k==i:
            continue
        if k%4==0:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta4,phi4)/(PD(i,k,L)**alpha)
        elif k%4==1:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta4,phi4)/(PD(i,k,L)**alpha)
        elif k%4==2:
            sumans += innerproduct_2vector_sphere(theta3,phi3,theta4,phi4)/(PD(i,k,L)**alpha)
        else:
            sumans += 1/(PD(i,k,L)**alpha)
    return sumans/4.0*2.0

def minimizationstate4optimization1D():
    L=int(sys.argv[1])
    initguess=np.random.random([6])*2*np.pi
    print("The initguess is {}".format(initguess))
    alpha=float(sys.argv[2])
    print(L,alpha)
    def state4optimization1Dwrapper(thetas):
        return state4optimization1D(thetas,L,alpha)
    res=minimize(state4optimization1Dwrapper,initguess,method='BFGS',)
    print("最小化的角度是{}".format(res.x))
    print("最小化的值是{}".format(res.fun))
    print(res)


def state5optimization1D(thetas,L,alpha):
    theta1=0
    phi1=0
    theta2,phi2,theta3,phi3,theta4,phi4,theta5,phi5=thetas
    sumans=0
    i=0
    for k in range(1,L):
        if k%5==0:
            sumans += 1/(PD(i,k,L)**alpha)
        elif k%5==1:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta2,phi2)/(PD(i,k,L)**alpha)
        elif k%5==2:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta3,phi3)/(PD(i,k,L)**alpha)
        elif k%5==3:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta4,phi4)/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta5,phi5)/(PD(i,k,L)**alpha)
    i=1
    for k in range(L):
        if k==i:
            continue
        if k%5==0:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta2,phi2)/(PD(i,k,L)**alpha)
        elif k%5==1:
            sumans += 1/(PD(i,k,L)**alpha)
        elif k%5==2:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta3,phi3)/(PD(i,k,L)**alpha)
        elif k%5==3:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta4,phi4)/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta5,phi5)/(PD(i,k,L)**alpha)
    i=2
    for k in range(L):
        if k==i:
            continue
        if k%5==0:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta3,phi3)/(PD(i,k,L)**alpha)
        elif k%5==1:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta3,phi3)/(PD(i,k,L)**alpha)
        elif k%5==2:
            sumans += 1/(PD(i,k,L)**alpha)
        elif k%5==3:
            sumans += innerproduct_2vector_sphere(theta3,phi3,theta4,phi4)/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta3,phi3,theta5,phi5)/(PD(i,k,L)**alpha)
    i=3
    for k in range(L):
        if k==i:
            continue
        if k%5==0:
            sumans += innerproduct_2vector_sphere(theta1,phi1,theta4,phi4)/(PD(i,k,L)**alpha)
        elif k%5==1:
            sumans += innerproduct_2vector_sphere(theta2,phi2,theta4,phi4)/(PD(i,k,L)**alpha)
        elif k%5==2:
            sumans += innerproduct_2vector_sphere(theta3,phi3,theta4,phi4)/(PD(i,k,L)**alpha)
        elif k%5==3:
            sumans += 1/(PD(i,k,L)**alpha)
        else:
            sumans += innerproduct_2vector_sphere(theta4,phi4,theta5,phi5)/(PD(i,k,L)**alpha)
    i=4
    



if __name__=='__main__':
   # different_alpha_correlators(10)
   minimizationsite2D1()


