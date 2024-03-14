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


def solve_classical_system(tend,tnum,h,y):
    t=np.linspace(0,tend,tnum)
    sol=solve_ivp(classical_eq,[0,tend],y,t_eval=t,args=(h,))

def theta_interpolated_function(sol):
    return interp1d(sol.t,sol.y[0],kind='cubic',fill_value='extrapolate')

def phi_interpolated_function(sol):
    return interp1d(sol.t,sol.y[1],kind='cubic',fill_value='extrapolate')

def system_G(t,y,)



