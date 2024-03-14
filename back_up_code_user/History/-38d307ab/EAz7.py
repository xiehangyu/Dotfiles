from trial import *

def drawing_entanglemetn(N,alpha,J,tend):
    H=simple_longrange_Hamiltonian(N,alpha,J)
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

if __name__=='__main__':
    for i in range(0,21):
        h=i/20.0
        drawing_entanglemetn(10,0.1,h,120)
    
