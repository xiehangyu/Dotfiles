from trial import *

def drawing_entanglemetn(N,alpha,J,tend):
    H=simple_longrange_Hamiltonian(N,alpha,J)
    print(H.eigenenergies())
    psi0=tensor([basis(2,0)]*N)
    tlist=np.linspace(0,tend,100)
    result=mesolve(H,psi0,tlist,[],[])
    times=result.times
    states=result.states
    entanglement=[]
    for i in range(len(states)):
        entanglement.append(entropy_vn(states[i].ptrace([j for j in range(N//2)])))
    plt.plot(times,entanglement)
    plt.show()

def drawing_in_one_picture(N,J,tend):
    ts=np.linspace(0,tend,4000)
    entropyss=[]
    for alpha in {0,0.1}:
        H=simple_longrange_HamiltonianOBC(N,alpha,J)
        psi0=tensor([basis(2,0)]*N)
        result=mesolve(H,psi0,ts,[],[])
        times=result.times
        states=result.states
        entanglements=[]
        for i in range(len(states)):
            entanglements.append(entropy_vn(states[i].ptrace([j for j in range(N//2)])))
        entropyss.append(entanglements)
        plt.plot(times,entanglements,label='OBC{}'.format(alpha))
        H=simple_longrange_Hamiltonian(N,alpha,J)
        psi0=tensor([basis(2,0)]*N)
        result=mesolve(H,psi0,ts,[],[])
        times=result.times
        states=result.states
        entanglements=[]
        for i in range(len(states)):
            entanglements.append(entropy_vn(states[i].ptrace([j for j in range(N//2)])))
        entropyss.append(entanglements)
        plt.plot(times,entanglements,label='PBC{}'.format(alpha))
    H=NearestNeighbor(N,J)
    psi0=tensor([basis(2,0)]*N)
    result=mesolve(H,psi0,ts,[],[])
    times=result.times
    states=result.states
    entanglements=[]
    for i in range(len(states)):
        entanglements.append(entropy_vn(states[i].ptrace([j for j in range(N//2)])))
    entropyss.append(entanglements)
    plt.plot(times,entanglements,label='NNPBC')
    H=NearestNeighborOBC(N,J)
    psi0=tensor([basis(2,0)]*N)
    result=mesolve(H,psi0,ts,[],[])
    times=result.times
    states=result.states
    entanglements=[]
    for i in range(len(states)):
        entanglements.append(entropy_vn(states[i].ptrace([j for j in range(N//2)])))
    entropyss.append(entanglements)
    plt.plot(times,entanglements,label='NNOBC')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('Entanglement')
    plt.show()
        

def newnormfactor(N,alpha):
    normfactor=0
    for i in range(1,N):
        normfactor+=(1+distance_PC(0,i,N))**(alpha)
    return normfactor
def simple_longrange_Hamiltonian(N,alpha,J=0.0):
    Z=sigmaz()/2
    X=sigmax()/2
    return_H=Qobj(dims=[[2]*N,[2]*N])
    for i in range(N):
        for j in range(i+1,N):
            return_H += tensor([identity(2)]*i+[Z]+[identity(2)]*(j-i-1)+[Z]+[identity(2)]*(N-j-1))/(1+distance_PC(i,j,N))**alpha
    for i in range(N):
        return_H += J*tensor([identity(2)]*i+[X]+[identity(2)]*(N-i-1))
    return return_H

def dynamical_phase_transition(h):
    N=12
    H2=NearestNeighbor(N,h)
    H=simple_longrange_Hamiltonian(N,0.1,h)
    psi0=tensor([basis(2,0)]*N)
    tlist=np.linspace(0,1000,1000)
    Z=sigmaz()
    Z=tensor([Z]+[qeye(2)]*(N-1))
    result=mesolve(H,psi0,tlist,[],[Z])
    result2=mesolve(H2,psi0,tlist,[],[Z])
    times=result.times
    expectation=result.expect[0]
    expectation2=result2.expect[0]
    plt.plot(times,expectation,times,expectation2)
    plt.legend(["Long Range","Nearest Neighbor"])
    plt.show()


if __name__=='__main__':
    for i in range(0,21):
        h=i/20.0
        drawing_entanglemetn(10,0.1,h,120)
    
