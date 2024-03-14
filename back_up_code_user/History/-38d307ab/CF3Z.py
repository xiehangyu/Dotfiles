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
    ts=np.linspace(0,tend,1000)
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
        

if __name__=='__main__':
    for i in range(0,21):
        h=i/20.0
        drawing_entanglemetn(10,0.1,h,120)
    
