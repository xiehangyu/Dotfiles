import numpy as np 
import matplotlib.pyplot as plt
import math
import scipy.special as ssc
from qutip import *
from scipy.optimize import minimize 

def rho(NA,M,t,N=1,J=1):
    #//M is the subsystem size, NA is the size of the environment
    #print("M=",M,"t=",t)
    rho=np.zeros((M+1,M+1))
    for i in range(M+1):
        for j in range(M+1):
            rho[i,j]=(np.cos(2*J*t*(i-j)/N))**(NA)*np.sqrt(ssc.comb(M,i)*ssc.comb(M,j))/2.0**M*np.exp(-((i-M/2.0)**2-(j-M/2.0)**2)*2*J*t/N)
    print(np.trace(rho))
    return rho


def entropy(rho):
    eig=np.linalg.eig(rho)
    eig=eig[0]
    eig=np.real(eig)
    eig=np.sort(eig)
    eig=eig[::-1]
    eig=eig[eig>0]
    return -np.sum(eig*np.log2(eig))

def plotting():
    xs=[]
    ys=[]
    for m in range(10,200,10):
        ys.append([])
        xs.append(str(m))
        for t in np.linspace(0,1E-5,500):
            ys[-1].append(entropy(rho(m,2,t,m+2)))
    for i in range(len(xs)):
        plt.plot(np.linspace(0,1E-5,500),ys[i])
        print("the rate is ",(ys[i][-1]-ys[i][0])/1E-5)
    plt.legend(xs)
    plt.show()

def plotting2():
    xs=[]
    ys=[]
    t=1E-4
    for m in range(2,90,1):
        xs.append(m)
        ys.append(entropy(rho(2*m,2,t,2*m+2))/t)
    fp=open("./plot_data/fitting.txt","w")
    for i in range(len(xs)):
        fp.write(str(xs[i])+" "+str(ys[i])+"\n")
    fp.close()
    xs=np.log(xs)
    ys=np.log(ys)
    plt.plot(xs,ys)
    print("the rate is ",(ys[-1]-ys[0])/(xs[-1]-xs[0]))
    plt.show()

def plotting5():
    xs=[]
    ys=[]
    t=1E-1
    tspace=50
    for m in range(10,110,10):
        print("m=",m)
        xs.append(str(m))
        ys.append([])
        for t in np.linspace(0,t,tspace):
            ys[-1].append(entropy(rho(200-m,m,t)))
    for i in range(len(xs)):
        plt.plot(np.linspace(0,t,tspace),ys[i])
        print("the rate is ",(ys[i][-1]-ys[i][0])/t)
    plt.legend(xs)
    plt.show()
def plotting4():
    xs=[]
    ys=[]
    t=1E-5
    for m in range(2,40,1):
        xs.append(m)
        ys.append(entropy2XXZ(t,2*m,m)/t)
    xs=np.log(xs)
    ys=np.log(ys)
    plt.plot(xs,ys)
    print("the rate is ",(ys[-1]-ys[0])/(xs[-1]-xs[0]))
    plt.show()

def writingfile():
    t=0.03
    fp=open("./plot_data/entropypermutation.txt","w")
    for m in range(10,210,10):
        entropyans=entropy(rho(400-m,m,t))
        fp.write(str(m)+" "+str(entropyans)+"\n")
    fp.close()



def XXZ_Hamiltonian(N,Delta=0.3):
    Njx=jmat(N/2,'x')
    Njy=jmat(N/2,'y')
    Njz=jmat(N/2,'z')
    H=Njx * Njx+Njy * Njy+Delta * Njz * Njz
    return H

def initialstate(N):
    state=Qobj(dims=[[N+1],[1]])
    for i in range(0,N+1):
        state += basis(N+1,i)*np.sqrt(ssc.comb(N,i)/2**N)
    #print(state.norm())
    state=state.unit()
    return state


def split_state(state,N,NA):
    NB=N-NA
    returnstate=Qobj(dims=[[NA+1,NB+1],[1,1]])
    for i in range(0,N+1):
        for j in range(max(0,i-NB),min(NA,i)+1):
            coef=np.sqrt(ssc.comb(NA,j)*ssc.comb(NB,i-j)/ssc.comb(N,i))
            returnstate += tensor(basis(NA+1,j),basis(NB+1,i-j))*coef*state[i,0]
    #print(returnstate.norm())
    returnstate=returnstate.unit()
    return returnstate

def entropy2XXZ(t,N,NA,Delta=0.3):
    H=XXZ_Hamiltonian(N,Delta)
    initial=initialstate(N)
    state=(-1j*t*H).expm()*initial
    state=state.unit()
    splitstate=split_state(state,N,NA)
    NAdensity=splitstate.ptrace(0)
    return entropy_vn(NAdensity,base=2)

def plotting3(t=100,t_total=100,Nmin=2,Nmax=30,Nstep=2):
    xs=[]
    ys=[]
    ts=np.linspace(0,t,t_total)
    for m in range(Nmin,Nmax,Nstep):
        xs.append(str(m))
        ys.append([])
        for t in ts:
            ys[-1].append(entropy2XXZ(t,2*m,m))
    for i in range(len(xs)):
        plt.plot(ts,ys[i])
        print("the rate is ",(ys[i][-1]-ys[i][0])/t)
    plt.legend(xs)
    plt.show()


def alltoallX(N):
    X=sigmax()/2
    returnops=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        returnops += tensor([identity(2)]*i+[X]+[identity(2)]*(N-i-1))
    return returnops

def alltoallY(N):
    Y=sigmay()/2
    returnops=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        returnops += tensor([identity(2)]*i+[Y]+[identity(2)]*(N-i-1))
    return returnops

def alltoallZ(N):
    Z=sigmaz()/2
    returnops=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        returnops += tensor([identity(2)]*i+[Z]+[identity(2)]*(N-i-1))
    return returnops

def alltoallplus(N):
    plus=sigmap()
    returnops=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        returnops += tensor([identity(2)]*i+[plus]+[identity(2)]*(N-i-1))
    return returnops

def totalspinop(N):
    X=alltoallX(N)
    Y=alltoallY(N)
    Z=alltoallZ(N)
    return X*X+Y*Y+Z*Z

def fourier_X(N,m):
    k=2*np.pi*m/N
    X=sigmax()/2
    returnops=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        returnops += tensor([identity(2)]*i+[X]+[identity(2)]*(N-i-1))*np.exp(1j*k*i)
    return returnops/N

def distance_PC(i,j,N):
    return min(abs(i-j),N-abs(i-j))

def NearestNeighbor(N):
    X=sigmax()
    returnop=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N-1):
        returnop += tensor([identity(2)]*i+[X]+[X]+[identity(2)]*(N-i-2))
    returnop += tensor([X]+[identity(2)]*(N-2)+[X])
    return returnop


def NearestNeighbor2(N):
    X=sigmax()
    Y=sigmay()
    Z=sigmaz()
    delta=0.7
    J=0.5
    returnop=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N-1):
        returnop += tensor([identity(2)]*i+[X]+[X]+[identity(2)]*(N-i-2))
        returnop += (1-delta)*tensor([identity(2)]*i+[Y]+[Y]+[identity(2)]*(N-i-2))
        returnop += tensor([identity(2)]*i+[Z]+[Z]+[identity(2)]*(N-i-2))
    returnop += tensor([X]+[identity(2)]*(N-2)+[X])
    returnop += tensor([Y]+[identity(2)]*(N-2)+[Y])
    returnop += (1-delta)*tensor([Z]+[identity(2)]*(N-2)+[Z])
    for i in range(N):
        returnop += J*tensor([identity(2)]*i+[Z]+[identity(2)]*(N-i-1))
    return returnop

def longrangeH2(N,alpha):
    Y=sigmay()/2
    Z=sigmaz()/2
    X=sigmax()/2
    delta=0.7
    J=0.5
    return_H=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        for j in range(i+1,N):
            return_H += tensor([identity(2)]*i+[X]+[identity(2)]*(j-i-1)+[X]+[identity(2)]*(N-j-1))/(1+distance_PC(i,j,N))**alpha
            return_H += (1-delta)*tensor([identity(2)]*i+[Y]+[identity(2)]*(j-i-1)+[Y]+[identity(2)]*(N-j-1))/(1+distance_PC(i,j,N))**alpha
            return_H += tensor([identity(2)]*i+[Z]+[identity(2)]*(j-i-1)+[Z]+[identity(2)]*(N-j-1))/(1+distance_PC(i,j,N))**alpha
        return_H += J*tensor([identity(2)]*i+[Z]+[identity(2)]*(N-i-1))
    return return_H/N**(1-alpha)

def simple_longrange_Hamiltonian(N,alpha):
    Y=sigmay()/2
    Z=sigmaz()/2
    X=sigmax()/2
    return_H=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        for j in range(i+1,N):
            return_H += tensor([identity(2)]*i+[X]+[identity(2)]*(j-i-1)+[X]+[identity(2)]*(N-j-1))/(1+distance_PC(i,j,N))**alpha
    return return_H/N**(1-alpha)

def evolution_simple_longrange(N,alpha,t):
    H=simple_longrange_Hamiltonian(N,alpha)
    initial=tensor([basis(2,0)]*N)
    expect_op=totalspinop(N)
    results=mesolve(H,initial,np.linspace(0,t,100),[],[])
    #plot the evolution of the expectation value
    times=results.times
    states=results.states
    expect_values=[]
    for state in states:
        expect_values.append(expect(expect_op,state))
    plt.plot(times,np.array(expect_values)/(N**2/4+N/2))
    plt.show()
    #plot the half-cut entanglement entropy of the system
    halfcut_entropies=[]
    for state in states:
        halfcut_entropies.append(entropy_vn(state.ptrace([i for i in range(N//2)]),base=2))
    reference_H=simple_longrange_Hamiltonian(N,0)
    reference_results=mesolve(reference_H,initial,np.linspace(0,t,100),[],[])
    reference_states=reference_results.states
    reference_halfcut_entropies=[]
    for state in reference_states:
        reference_halfcut_entropies.append(entropy_vn(state.ptrace([i for i in range(N//2)]),base=2))
    plt.plot(times,halfcut_entropies)
    plt.plot(times,np.array(halfcut_entropies)-np.array(reference_halfcut_entropies))
    plt.show()

def drawing_collapse(alpha,t):
    for N in range(5,14):
        H=simple_longrange_Hamiltonian(N,alpha)
        initial=tensor([basis(2,0)]*N)
        expect_op=totalspinop(N)
        results=mesolve(H,initial,np.linspace(0,t,100),[],[])
        #plot the evolution of the expectation value
        times=results.times
        states=results.states
        expect_values=[]
        for state in states:
            expect_values.append(expect(expect_op,state))
        plt.plot(times,(expect_values[0]-np.array(expect_values))/(N**2))
    plt.legend([str(i) for i in range(5,14)])
    plt.show()

def simple_longrange_Hamiltonian2(N,alpha):
    Y=sigmay()/2
    Z=sigmaz()/2
    X=sigmax()/2
    return_H=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        for j in range(i+1,N):
            return_H += tensor([identity(2)]*i+[X]+[identity(2)]*(j-i-1)+[X]+[identity(2)]*(N-j-1))/distance_PC(i,j,N)**alpha
    return return_H/N**(1-alpha)

def evolution_simple_longrange2(N,alpha,t,ifdrawing=True):
    H=simple_longrange_Hamiltonian(N,alpha)
    initial=tensor([basis(2,0)]*N)
    expect_op=totalspinop(N)
    results=mesolve(H,initial,np.linspace(0,t,100),[],[])
    #plot the evolution of the expectation value
    times=results.times
    states=results.states
    expect_values=[]
    for state in states:
        expect_values.append(expect(expect_op,state))
    plt.plot(times,np.array(expect_values)/(N**2/4+N/2))
    if ifdrawing:
        plt.show()
    #plot the half-cut entanglement entropy of the system
    halfcut_entropies=[]
    for state in states:
        halfcut_entropies.append(entropy_vn(state.ptrace([i for i in range(N//2)]),base=2))
    fp=open("./plot_data/HalfcutEntropy_N_{}_alpha_{}.txt".format(N,alpha),'w')
    for i in range(len(times)):
        fp.write(str(times[i])+" "+str(halfcut_entropies[i])+"\n")
    reference_H=hkcomponent(N,0,alpha)
    reference_results=mesolve(reference_H,initial,np.linspace(0,t,100),[],[])
    reference_states=reference_results.states
    reference_halfcut_entropies=[]
    for state in reference_states:
        reference_halfcut_entropies.append(entropy_vn(state.ptrace([i for i in range(N//2)]),base=2))
    plt.plot(times,halfcut_entropies)
    plt.plot(times,np.array(halfcut_entropies)-np.array(reference_halfcut_entropies))
    if ifdrawing:
        plt.show()

def evolution_simple_longrange2nearest(N,alpha,t,ifdrawing=True):
    H=NearestNeighbor(N)
    initial=tensor([basis(2,0)]*N)
    expect_op=totalspinop(N)
    results=mesolve(H,initial,np.linspace(0,t,100),[],[])
    #plot the evolution of the expectation value
    times=results.times
    states=results.states
    expect_values=[]
    for state in states:
        expect_values.append(expect(expect_op,state))
    plt.plot(times,np.array(expect_values)/(N**2/4+N/2))
    if ifdrawing:
        plt.show()
    #plot the half-cut entanglement entropy of the system
    halfcut_entropies=[]
    for state in states:
        halfcut_entropies.append(entropy_vn(state.ptrace([i for i in range(N//2)]),base=2))
    fp=open("./plot_data/HalfcutEntropyNearest_N_{}.txt".format(N),'w')
    for i in range(len(times)):
        fp.write(str(times[i])+" "+str(halfcut_entropies[i])+"\n")
    reference_H=hkcomponent(N,0,alpha)
    reference_results=mesolve(reference_H,initial,np.linspace(0,t,100),[],[])
    reference_states=reference_results.states
    reference_halfcut_entropies=[]
    for state in reference_states:
        reference_halfcut_entropies.append(entropy_vn(state.ptrace([i for i in range(N//2)]),base=2))
    plt.plot(times,halfcut_entropies)
    plt.plot(times,np.array(halfcut_entropies)-np.array(reference_halfcut_entropies))
    if ifdrawing:
        plt.show()
def hkcomponent(N,m,alpha):
    k=2*np.pi*m/N
    returnops=fourier_X(N,m)*fourier_X(N,-m)*N**alpha
    prefactor=0
    for i in range(N):
        prefactor += np.exp(-1j*k*i)/(1+distance_PC(i,0,N)**alpha)
    return returnops*prefactor

def construct_total_symmetry_state(N,m):
    states=tensor([basis(2,1)]*N)
    plus_op=alltoallplus(N)
    for i in range(m):
        states=plus_op*states 
    states=states.unit()
    return states

def func_need_optimal(x):
    N=(len(x)-1)*2
    state=Qobj(dims=[[N+1],[1]])
    for i in range(N//2):
        state+=x[i]*basis(N+1,i)+x[i]*basis(N+1,N-i)
    state += x[N//2]*basis(N+1,N//2)
    state=state.unit()
    splitstate=split_state(state,N,N//2)
    return -entropy_vn(splitstate.ptrace([0]),base=2)

def optimize_arg(N):
    x0=np.ones(N//2+1)
    res=minimize(func_need_optimal,x0,method='nelder-mead',options={'xtol':1e-8,'disp':True})
    print("Maximum entropy is "+str(-res.fun))
    print("The optimal state is "+str(res.x))
    fp=open("optimal_state{}.txt".format(N),'w')
    fp.write(str(res.x))
    fp.close()



def random_initial_state(N):
    states=Qobj(dims=[[2]*N,[1]*N])
    for i in range(N+1):
        states+=np.random.rand()*construct_total_symmetry_state(N,i)
    states=states.unit()
    print("The expectation value of total spin is "+str(expect(alltoallZ(N),states)))
    return states.unit()

def upperbound(N,Squar,NA):
    S_o=N/2
    deltaS=S_o*(S_o+1)-Squar
    #lamb=deltaS/(2*S_o)
    lamb=deltaS/(S_o*(S_o+1))
    np.seterr(all='raise')
    try:
        lamb_sq=np.sqrt(lamb)
    except:
        print("error")
        print(lamb)
        lamb=0
        lamb_sq=np.sqrt(lamb)
    entropy_origin=np.log2(NA+1)
    try:
        entropy_diff=lamb_sq*np.log2(2**NA-1)-lamb_sq*np.log2(lamb_sq)-(1-lamb_sq)*np.log2(1-lamb_sq)
        #entropy_diff=lamb_sq*np.log2(NA)-lamb_sq*np.log2(lamb_sq)-(1-lamb_sq)*np.log2(1-lamb_sq)
    except:
        print(lamb_sq)
        entropy_diff=lamb_sq*np.log2(2**NA-1)
        #entropy_diff=lamb_sq*np.log2(NA)
    return entropy_origin+entropy_diff

def Fobbnes_inequality(NA,T):
    d=2**NA
    if T>1-1/d:
        return np.log2(d)
    else:
        if T==0:
            return 0
        return T*np.log2(d-1)-T*np.log2(T)-(1-T)*np.log2(1-T)

def read_state(N):
    fp=open("optimal_state{}.txt".format(N),'r')
    state=Qobj(dims=[[2]*N,[1]*N])
    x_vec=fp.read().split()
    x_vec=x_vec[1:N//2+2]
    x_vec[-1]=x_vec[-1][:-1]
    for i in range(N//2):
        state+=float(x_vec[i])*construct_total_symmetry_state(N,i)+float(x_vec[i])*construct_total_symmetry_state(N,N-i)
    state+=float(x_vec[N//2])*construct_total_symmetry_state(N,N//2)
    state=state.unit()
    return state.unit()

def state_project(N,state):
    return_state=Qobj(dims=[[2]*N,[1]*N])
    for i in range(N+1):
        temp_state=construct_total_symmetry_state(N,i)
        return_state=return_state+temp_state*state.overlap(temp_state)
    return return_state.unit() 


def project_operator(N):
    return_project=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N+1):
        state=construct_total_symmetry_state(N,i)
        return_project+=state*state.dag()
    return return_project

def sqrttracedist(rho1,rho2):
    rho1_sqrt=rho1.sqrtm()
    rho2_sqrt=rho2.sqrtm()
    return 2*tracedist(rho1_sqrt,rho2_sqrt)


def overlapstates_totalangularoperator(rho,N,alpha=0.1,if_drawing=True):
    op=totalspinop(N)
    states=op.eigenstates()
    fp=open("./plot_data/overlapstates_N_{}_alpha_{}.txt".format(N,alpha),'w')
    overlaps_c=[]
    totalindex=len(states[0])
    for i in range(len(states[0])):
        overlapvalue=abs(rho.overlap(states[1][totalindex-i-1]))**2
        fp.write("{}  {}  {}\n".format(i,states[0][totalindex-i-1],overlapvalue))
        overlaps_c.append(overlapvalue)
    fp.close()
    energies_c=states[0][::-1]
    energy_count=[]
    count_c=[]
    previous=energies_c[0]
    energy_count.append(previous)
    temp_count=[]
    for i in range(len(overlaps_c)):
        if np.abs(energies_c[i]-previous)>1E-3:
            count_c.append(np.sum(np.array(temp_count)))
            temp_count=[]
            temp_count.append(overlaps_c[i])
            energy_count.append(energies_c[i])
        else:
            temp_count.append(overlaps_c[i])
        previous=energies_c[i]
    count_c.append(np.sum(np.array(temp_count)))
    fp=open("plot_data/overlapstatescounting_N_{}_alpha_{}.txt".format(N,alpha),'w')
    for i in range(len(count_c)):
        fp.write("{}  {}  {}\n".format(i,energy_count[i],count_c[i]))
    fp.close()
    if if_drawing:
        print(len(count_c))
        print(np.sum(np.array(count_c)))
        print(np.sum(np.array(overlaps_c)))
        plt.bar(range(len(count_c)),count_c)
        plt.show()

def largesteigenvalues(rho,N,alpha=0.1,if_drawing=True):
    op=project_operator(N//2)
    temprho=rho.ptrace([i for i in range(N//2)])
    eigenvalues=temprho.eigenstates()
    eigenstates_c=eigenvalues[1]
    eigenvalues=eigenvalues[0]
    eigenvalues=eigenvalues[::-1]
    eigenstates_c=eigenstates_c[::-1]
    fp=open("./plot_data/largesteigenvalues_N_{}_alpha_{}.txt".format(N,alpha),'w')
    drawing_eigs=[]
    drawing_overlaps=[]
    for i in range(np.min([50,2**(N//2)])):
        drawing_overlap=expect(op,eigenstates_c[i])
        fp.write("{}  {}  {}\n".format(i,eigenvalues[i],drawing_overlap))
        drawing_eigs.append(eigenvalues[i])
        drawing_overlaps.append(drawing_overlap)
    fp.close()
    if if_drawing:
        plt.plot(np.array(drawing_eigs))
        plt.show()
        plt.plot(np.array(drawing_overlaps))
        plt.show()
def state_in_spinS(N,S):
    #Here S is the total spin, 0 \leq S \leq N/2
    beginindex=N/2-S
    if beginindex==0:
        repeat=1
    else:
        repeat=ssc.comb(N,beginindex)-ssc.comb(N,beginindex-1)
    return repeat*(N+1-2*beginindex)

def rank_of_states(results,cutoff=1E-7):
    times=results.times
    states=results.states
    ranks=[]
    for i in range(len(times)):
        state=states[i]
        state=state.ptrace([i for i in range(len(state.dims[0]))//2])
        eigenvalues=state.eigenenergies()
        eigenvalues=eigenvalues[::-1]
        eigenvalues=eigenvalues/np.sum(eigenvalues) 
        rank=0
        for j in range(len(eigenvalues)):
            if eigenvalues[j]>cutoff:
                rank+=1
        ranks.append(rank)
    plt.plot(times,ranks)
    plt.show()

def trial_eigenvalues(N,alpha=0.1):
    #N=input("please input the size of the system\n")
    #N=int(N)
    H=longrangeH2(N,alpha)
    if alpha==-1:
        H=NearestNeighbor2(N)
    opts = Options(nsteps=10000)
    initial=tensor([basis(2,0)]*N)
    results=mesolve(H,initial,np.linspace(0,200,100),[],[],options=opts)
    rho=results.states[-1]
    largesteigenvalues(rho,N,alpha,True)
    overlapstates_totalangularoperator(rho,N,alpha,True)

def sqrthilbertdist(rho1,rho2):
    rho1_sqrt=rho1.sqrtm()
    rho2_sqrt=rho2.sqrtm()
    return np.sqrt(hilbert_dist(rho1_sqrt,rho2_sqrt))

def evolution_longrange_tracedistance(N,alpha,t,if_drawing=True):
    #H=simple_longrange_Hamiltonian(N,alpha)
    H=longrangeH2(N,alpha)
    #P=project_operator(N)
    initial=tensor([basis(2,1)]*N)
    #initial=random_initial_state(N)
    results=mesolve(H,initial,np.linspace(0,t,400),[],[])
    times=results.times
    states=results.states
    tracedistances=[]
    sqrttracedistances=[]
    sqrthilbertdistances=[]
    for state in states:
        #referencestate=P*state
        referencestate=state_project(N,state)
        referencestate=referencestate.unit()
        tracedistances.append(2*tracedist(state.ptrace([i for i in range(N//2)]),referencestate.ptrace([i for i in range(N//2)])))
        sqrttracedistances.append(sqrttracedist(state.ptrace([i for i in range(N//2)]),referencestate.ptrace([i for i in range(N//2)])))
        sqrthilbertdistances.append(sqrthilbertdist(state.ptrace([i for i in range(N//2)]),referencestate.ptrace([i for i in range(N//2)])))
    plt.rcParams.update({'font.size': 25})
    plt.plot(times,tracedistances)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$D(\rho,\rho_{\mathrm{sym}})$")
    if if_drawing:
        plt.show()
    plt.plot(times,np.exp(np.array(tracedistances)))
    plt.xlabel(r"$t$")
    plt.ylabel(r"$Exp(D(\rho,\rho_{\mathrm{sym}}))$")
    if if_drawing:
        plt.show()
    entropy_halfcut=[]
    entropy_up=[]
    NA=N//2
    for i in range(len(times)):
        entropy_halfcut.append(entropy_vn(states[i].ptrace([i for i in range(NA)]),base=2))
        entropy_up.append(np.log2(NA+1)+Fobbnes_inequality(NA,tracedistances[i]))
    plt.plot(times,entropy_halfcut)
    plt.plot(times,entropy_up)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$S$")
    plt.legend([r"$S_{\mathrm{halfcut}}$",r"$S_{\mathrm{upperbound}}$"])
    if if_drawing:
        plt.show()
    plt.plot(times,sqrttracedistances)
    plt.plot(times,sqrthilbertdistances)
    plt.plot(times,tracedistances)
    fp=open("./plot_data/distancelongrange_N_{}_alpha_{}.txt".format(N,alpha),'w')
    for i in range(len(times)):
        fp.write(str(times[i])+" "+str(sqrttracedistances[i])+" "+str(sqrthilbertdistances[i])+" "+str(tracedistances[i])+"\n")
    fp.close()
    plt.xlabel(r"$t$")
    plt.ylabel(r"sqrt distance")
    plt.legend([r"$\sqrt{D(\rho,\rho_{\mathrm{sym}})}$",r"$\sqrt{D_{\mathrm{Hilbert}}(\rho,\rho_{\mathrm{sym}})}$","tracedistance"])
    if if_drawing:
        plt.show()
    distratios=[]
    for i in range(len(sqrthilbertdistances)):
        if sqrthilbertdistances[i]==0:
            distratios.append(0)
        else:
            distratios.append(sqrttracedistances[i]/sqrthilbertdistances[i])
    plt.plot(times,distratios)
    plt.xlabel(r"$t$")
    plt.ylabel("ratios")
    if if_drawing:
        plt.show()
    expect_op=totalspinop(N)
    expect_values=[]
    for state in states:
        expect_values.append(expect(expect_op,state))
    lambdas=[]
    lambdas2=[]
    expop=[]
    S_o=N/2
    for i in range(len(times)):
        lambdas.append(np.sqrt((S_o*(S_o+1)-expect_values[i])/(2*S_o)))
        lambdas2.append(np.sqrt((S_o*(S_o+1)-expect_values[i])/(S_o*(S_o+1))))
        expop.append(S_o*(S_o+1)-expect_values[i])
    plt.plot(times,expop)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\langle S^2\rangle$")
    if if_drawing:
        plt.show()
    plt.plot(times,tracedistances,times,lambdas,times,lambdas2)
    plt.legend([r"$D(\rho,\rho_{\mathrm{sym}})$",r"$\lambda$",r"$\lambda^2$"])
    plt.xlabel(r"$t$")
    if if_drawing:
        plt.show()

def evolution_longrange_tracedistanceNearestneighbor(N,alpha,t,if_drawing=True):
    #H=NearestNeighbor(N)
    H=NearestNeighbor2(N)
    #P=project_operator(N)
    initial=tensor([basis(2,1)]*N)
    #initial=random_initial_state(N)
    results=mesolve(H,initial,np.linspace(0,t,400),[],[])
    times=results.times
    states=results.states
    tracedistances=[]
    sqrttracedistances=[]
    sqrthilbertdistances=[]
    for state in states:
        #referencestate=P*state
        referencestate=state_project(N,state)
        referencestate=referencestate.unit()
        tracedistances.append(2*tracedist(state.ptrace([i for i in range(N//2)]),referencestate.ptrace([i for i in range(N//2)])))
        sqrttracedistances.append(sqrttracedist(state.ptrace([i for i in range(N//2)]),referencestate.ptrace([i for i in range(N//2)])))
        sqrthilbertdistances.append(sqrthilbertdist(state.ptrace([i for i in range(N//2)]),referencestate.ptrace([i for i in range(N//2)])))
    plt.rcParams.update({'font.size': 25})
    plt.plot(times,tracedistances)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$D(\rho,\rho_{\mathrm{sym}})$")
    if if_drawing:
        plt.show()
    plt.plot(times,np.exp(np.array(tracedistances)))
    plt.xlabel(r"$t$")
    plt.ylabel(r"$Exp(D(\rho,\rho_{\mathrm{sym}}))$")
    if if_drawing:
        plt.show()
    entropy_halfcut=[]
    entropy_up=[]
    NA=N//2
    for i in range(len(times)):
        entropy_halfcut.append(entropy_vn(states[i].ptrace([i for i in range(NA)]),base=2))
        entropy_up.append(np.log2(NA+1)+Fobbnes_inequality(NA,tracedistances[i]))
    plt.plot(times,entropy_halfcut)
    plt.plot(times,entropy_up)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$S$")
    plt.legend([r"$S_{\mathrm{halfcut}}$",r"$S_{\mathrm{upperbound}}$"])
    if if_drawing:
        plt.show()
    fp=open("./plot_data/tracedistance.txt",'w')
    for i in range(len(times)):
        fp.write(str(times[i])+" "+str(tracedistances[i])+"\n")
    fp.close()
    expect_op=totalspinop(N)
    expect_values=[]
    for state in states:
        expect_values.append(expect(expect_op,state))
    lambdas=[]
    lambdas2=[]
    expop=[]
    S_o=N/2
    for i in range(len(times)):
        lambdas.append(np.sqrt((S_o*(S_o+1)-expect_values[i])/(2*S_o)))
        lambdas2.append(np.sqrt((S_o*(S_o+1)-expect_values[i])/(S_o*(S_o+1))))
        expop.append(S_o*(S_o+1)-expect_values[i])
    plt.plot(times,expop)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\langle S^2\rangle$")
    if if_drawing:
        plt.show()
    plt.plot(times,tracedistances,times,lambdas,times,lambdas2)
    plt.legend([r"$D(\rho,\rho_{\mathrm{sym}})$",r"$\lambda$",r"$\lambda^2$"])
    plt.xlabel(r"$t$")
    if if_drawing:
        plt.show()
    plt.plot(times,sqrttracedistances)
    plt.plot(times,sqrthilbertdistances)
    plt.plot(times,tracedistances)
    fp=open("./data/distancenearest_N_{}.txt".format(N),'w')
    for i in range(len(times)):
        fp.write(str(times[i])+" "+str(sqrttracedistances[i])+" "+str(sqrthilbertdistances[i])+" "+str(tracedistances[i])+"\n") 
    fp.close()
    plt.xlabel(r"$t$")
    plt.ylabel(r"sqrt distance")
    plt.legend([r"$\sqrt{D(\rho,\rho_{\mathrm{sym}})}$",r"$\sqrt{D_{\mathrm{Hilbert}}(\rho,\rho_{\mathrm{sym}})}$","tracedistance"])
    if if_drawing:
        plt.show()
    distratios=[]
    for i in range(len(sqrthilbertdistances)):
        if sqrthilbertdistances[i]==0:
            distratios.append(0)
        else:
            distratios.append(sqrttracedistances[i]/sqrthilbertdistances[i])
    plt.plot(times,distratios)
    plt.xlabel(r"$t$")
    plt.ylabel("ratios")
    if if_drawing:
        plt.show()

def compare_output(N,alpha,ts):
    H=simple_longrange_Hamiltonian(N,alpha)
    initial=tensor([basis(2,0)]*N)
    expect_op=totalspinop(N)
    for t in ts:
        #results=mesolve(H,initial,np.linspace(0,t,100),[],[])
        #state=results.states[-1]
        state=(-1j*H*t).expm()*initial
        print("{} {}".format(t,N/2*(N/2+1)-expect(expect_op,state)))
        print("{} {}".format(t,entropy_vn(state.ptrace([i for i in range(N//2)]),base=2)))

def evolution_longrange_simplebound(N,alpha,t):
    H=simple_longrange_Hamiltonian(N,alpha)
    initial=tensor([basis(2,1)]*N)
    #initial=random_initial_state(N)
    expect_op=totalspinop(N)
    results=mesolve(H,initial,np.linspace(0,t,100),[],[])
    #plot the evolution of the expectation value
    times=results.times
    states=results.states
    expect_values=[]
    for state in states:
        expect_values.append(expect(expect_op,state))
    plt.rcParams.update({'font.size': 25})
    plt.plot(times,np.array(expect_values)/(N**2/4+N/2))
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\frac{\langle S^2\rangle}{N/2(N/2+1)}$")
    plt.show()
    #plot the half-cut entanglement entropy of the system
    halfcut_entropies=[]
    for state in states:
        halfcut_entropies.append(entropy_vn(state.ptrace([i for i in range(N//2)]),base=2))
    reference_H=simple_longrange_Hamiltonian(N,0)
    reference_results=mesolve(reference_H,initial,np.linspace(0,t,len(times)),[],[])
    reference_states=reference_results.states
    reference_halfcut_entropies=[]
    for state in reference_states:
        reference_halfcut_entropies.append(entropy_vn(state.ptrace([i for i in range(N//2)]),base=2))
    plt.plot(times,halfcut_entropies)
    plt.show()
    plt.plot(times,np.array(halfcut_entropies)-np.array(reference_halfcut_entropies))
    plt.show()
    entropy_up=[]
    for expect_va in expect_values:
        entropy_up.append(upperbound(N,expect_va,N//2))
    plt.plot(times,halfcut_entropies,times,entropy_up)
    plt.legend(['entropy','upperbound'])
    plt.xlabel(r"$t$")
    plt.ylabel(r"$S_{\rho}$")
    plt.show()
    plt.plot(times,np.exp(np.array(halfcut_entropies)))
    plt.show()
    plt.plot(times,np.exp(np.array(entropy_up)))
    plt.show()
    fp=open("./data/entropy_up.txt",'w')
    for i in range(len(times)):
        fp.write(str(times[i])+' '+str(entropy_up[i])+'\n')
    fp=open("./data/totalspin.txt",'w')
    for i in range(len(times)):
        fp.write(str(times[i])+' '+str(expect_values[i])+'\n')






if __name__=='__main__':
    #evolution_longrange_tracedistance(12,0.1,200,False)
    #evolution_longrange_tracedistanceNearestneighbor(12,0.1,200,False)
    #evolution_simple_longrange2(20,0.1,200,False)
    #evolution_simple_longrange2nearest(20,0.1,200,False)
    for i in range(6,24,2):
        trial_eigenvalues(i,-1)
