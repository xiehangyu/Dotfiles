import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

def period_length(i,N):
    return np.min([i,N-i])

def f_k(N,m,alpha):
    normfactor=1.0/np.power(N,1-alpha)
    k=2*np.pi/N*m
    sum_result=0
    for i in range(1,N):
        sum_result += normfactor*np.exp(-1j*k*i)/np.power(period_length(i,N)+1,alpha)
    return sum_result

def Jks(N,alpha):
    fks=[]
    for m in range(N):
        fks.append(f_k(N,m,alpha))
    return np.array(fks)

def classical_eq(t,y,h,J):
    theta,phi=y
    dthetadt=2*J*np.sin(theta)*np.cos(phi)*np.sin(phi)
    dphidt=-2*h+2*J*np.cos(theta)*np.cos(phi)*np.cos(phi)
    return [dthetadt,dphidt]


def solve_classical_system(tend,tnum,h,J,y):
    t=np.linspace(0,tend,tnum)
    sol=solve_ivp(classical_eq,[0,tend],y,t_eval=t,args=(h,J))

def theta_interpolated_function(sol):
    return interp1d(sol.t,sol.y[0],kind='cubic',fill_value='extrapolate')

def phi_interpolated_function(sol):
    return interp1d(sol.t,sol.y[1],kind='cubic',fill_value='extrapolate')

def Fourier_transformation(N):
    return_U=np.zeros((N,N),dtype=complex)
    for i in range(N):
        for m in range(N):
            return_U[i,m]=np.exp(2*np.pi*1j*i*m/N)*np.sqrt(2/N)
    return_U=np.kron(return_U,np.eye(2))
    return return_U

def truncation_matrix(matrix,NA):
    return matrix[:2*NA,:2*NA]

def symplectic_matrix(N):
    return np.kron(np.eye(N),np.array([[0,1],[-1,0]]))

def calculate_entropy(matrix):
    matrix_new=1j*matrix
    eigenalues=np.linalg.eigvals(matrix_new)
    entropy=0
    for i in eigenvalues:
        if i>0:
            entropy += (i+1)/2*np.log2((i+1)/2)-(i-1)/2*np.log2((i-1)/2)
    return entropy

class evolution_of_functions:
    def __init__(self,alpha,N,theta,phi,tend,tnum,h):
        self.theta=theta
        self.phi=phi
        self.N=N
        self.tend=tend
        self.tnum=tnum
        self.h=h
        self.alpha=alpha
        self.Js=Jks(N,alpha)
        self.sol=solve_classical_system(tend,tnum,h,self.Js[0],[theta,phi])
        self.theta_interpolated_function=theta_interpolated_function(self.sol)
        self.phi_interpolated_function=phi_interpolated_function(self.sol)
        self.J=self.Js[0]
        self.U=Fourier_transformation(N)
    def system_G(self,t,y,Jk):
        J=self.J
        Gqq,Gpp,Gqp=y
        theta=self.theta_interpolated_function(t)
        phi=self.phi_interpolated_function(t)
        dGqqdt=4*Jk*np.cos(theta)*np.cos(phi)*np.sin(phi)*Gqq+4*(J*np.cos(phi)**2-Jk*np.sin(phi)**2)*Gqp
        dGppdt=-4*(J*np.cos(phi)**2-Jk*np.cos(theta)**2*np.cos(phi)**2)*Gqp-4*Jk*np.cos(theta)*np.cos(phi)*np.sin(phi)*Gpp
        dGqpdt=2*(J*np.cos(phi)**2-Jk*np.sin(phi)**2)*Gpp-2*(J*np.cos(phi)**2-Jk*np.cos(theta)**2*np.cos(phi)**2)*Gqq
        return [dGqqdt,dGppdt,dGqpdt]

    def solve_system_G(self,Jk):
        t=np.linspace(0,self.tend,self.tnum)
        y=[1/2,1/2,0]
        sol=solve_ivp(self.system_G,[0,self.tend],y,t_eval=t,args=(Jk,))
        return sol.y
    
    def evolution_of_G_kspace(self):
        Gs=[]
        for Jk in self.Js:
            Gs.append(self.solve_system_G(Jk))
        self.Kmatrixs=[]
        self.ts=np.linspace(0,self.tend,self.tnum)
        for i in range(self.tnum):
            new_matrix=np.zeros((2*self.N,2*self.N),dtype=complex)
            for m in range(self.N):
                new_matrix[2*m,2*m]=Gs[m][0,i]
                new_matrix[2*m+1,2*m+1]=Gs[m][1,i]
                new_matrix[2*m,2*m+1]=Gs[m][2,i]
                new_matrix[2*m+1,2*m]=Gs[m][2,i]
            self.Kmatrixs.append(new_matrix)
    def G_position_matrix(self):
        self.Pmatrixs=[]
        for i in range(self.tnum):
            self.Pmatrixs.append(np.dot(self.U,np.dot(self.Kmatrixs[i],self.U.T.conj())))
    def Jmatrix(self):
        symplectic=symplectic_matrix(self.N)
        self.Jmatrixs=[]
        for i in range(self.tnum):
            self.Jmatrixs.append(np.dot(-self.Pmatrixs[i],symplectic))
    def Subsystemmatrix(self):
        self.submatrixs=[]
        for i in range(self.tnum):
            self.submatrixs.append(truncation_matrix(self.Jmatrixs[i],self.N))
    
    def allinall(self):
        self.evolution_of_G_kspace()
        self.G_position_matrix()
    








