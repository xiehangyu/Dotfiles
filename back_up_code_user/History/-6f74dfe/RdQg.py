import numpy as np  
import matplotlib.pyplot as plt  
import scipy 
from ncon import ncon
from decomposition import kron_product


def diagonal_phi_loop(W,Nv):
    Wdouble=np.square(W)
    Wdouble=Wdouble.reshape([2,2,2,2])
    start=Wdouble
    for i in range(1,Nv):
        start=ncon([start,Wdouble],[[-j for j in range(1,i+1)]+[1]+[-j for j in range(i+4,2*i+4)]+[-2*i-4],[-i-1,-i-2,-i-3,1]])
    id=np.eye(2)
    return ncon([start,id],[[-j for j in range(1,Nv+1)]+[1]+[-(3*Nv+1-j) for j in range(Nv+1,2*Nv+1)]+[2],[1,2]]).reshape([2**Nv,2**Nv])

def offdiagonal_phi_loop(W,Nv):
    Wdouble=np.square(W)
    Wdouble=Wdouble.reshape([2,2,2,2])
    start=Wdouble
    for i in range(1,Nv):
        start=ncon([start,Wdouble],[[-j for j in range(1,i+1)]+[1]+[-j for j in range(i+4,2*i+4)]+[-2*i-4],[-i-1,-i-2,-i-3,1]])
    Z=np.array([[1,0],[0,-1]])
    return ncon([start,Z],[[-j for j in range(1,Nv+1)]+[1]+[-j for j in range(Nv+1,2*Nv+1)]+[2],[1,2]]).reshape([2**Nv,2**Nv])


def diagonal_symmetry(W,Nv):
    loop_transfer=diagonal_phi_loop(W,Nv)
    Z=np.array([[1,0],[0,-1]])
    id=np.eye(2)
    projector=(kron_product([id]*Nv)+kron_product([Z]*Nv))/2.0
    return projector@loop_transfer@projector.T

def offdiagonal_symmetry(W,Nv):
    loop_transfer=offdiagonal_phi_loop(W,Nv)
    Z=np.array([[1,0],[0,-1]])
    id=np.eye(2)
    projector=(kron_product([id]*Nv)+kron_product([Z]*Nv))/2.0
    return projector@loop_transfer@projector.T

def diagonal_asymmetry(W,Nv):
    loop_transfer=diagonal_phi_loop(W,Nv)
    Z=np.array([[1,0],[0,-1]])
    id=np.eye(2)
    projector=(kron_product([id]*Nv)-kron_product([Z]*Nv))/2.0
    return projector@loop_transfer@projector.T

def offdiagonal_asymmetry(W,Nv):
    loop_transfer=offdiagonal_phi_loop(W,Nv)
    Z=np.array([[1,0],[0,-1]])
    id=np.eye(2)
    projector=(kron_product([id]*Nv)-kron_product([Z]*Nv))/2.0
    return projector@loop_transfer@projector.T

def create_W_matrix1(coef):
    W=np.array([[1,0,0,coef**2],[0,coef**2,coef**2,0],[0,coef**2,coef**2,0],[coef**2,0,0,coef**4]])
    return W
    
def create_W_matrix2(coef):
    t1=np.cos(coef)
    t2=np.sin(coef)
    W=np.array([[t1,0,0,t2],[0,t2,t1,0],[0,t2,t1,0],[t1,0,0,t2]])
    return W

def create_W_matrix3(alpha,beta):
    W=np.array([[alpha,0,0,1-alpha],[0,beta,1-beta,0],[0,1-alpha,alpha,0],[1-beta,0,0,beta]])
    W=np.sqrt(W)
    return W

def path1(coef):
    return coef,coef

def theoritical_diagonal_firsttwo(alpha,beta,Nv):
    if alpha==0 or beta==0:
        return 0,0
    item1=(1+beta/alpha)**Nv*(1-alpha)*(1-beta)/(alpha*beta)+1-(1-alpha)*(1-beta)/(alpha*beta)+(beta/alpha)**Nv*(alpha+beta-1)/(alpha*beta)
    item1=item1*alpha**Nv
    item2=(1+beta/alpha)**(Nv-2)*(1-beta/alpha)**2*(1-alpha)*(1-beta)/(alpha*beta)+1-(1-alpha)*(1-beta)/(alpha*beta)+(beta/alpha)**Nv*(alpha+beta-1)/(alpha*beta)
    item2=item2*alpha**Nv
    return item1,item2

def certification_or_theoritical_diagonal(alpha,beta,Nv):
    W=create_W_matrix3(alpha,beta)
    transferloop=diagonal_symmetry(W,Nv)
    m1s=np.sort(np.abs(np.linalg.eigvals(transferloop)))[::-1]
    item1,item2=theoritical_diagonal_firsttwo(alpha,beta,Nv)
    refer1,refer2=m1s[0],m1s[1]
    print("Prediction:",item1,item2)
    print("Numerical:",refer1,refer2)

def phase_boundary_overlap(coefs):
    satisfying_point_x=[]
    satisfying_point_y=[]
    Nv=12
    for alpha,beta in coefs:
        W=create_W_matrix3(alpha,beta)
        m1=np.sort(np.abs(np.linalg.eigvals(diagonal_symmetry(W,Nv))))[-1]
        m3=np.sort(np.abs(np.linalg.eigvals(offdiagonal_symmetry(W,Nv))))[-1]
        if np.abs(1-m3/m1)<=2E-6:
            satisfying_point_x.append(alpha)
            satisfying_point_y.append(beta)
    plt.scatter(satisfying_point_x,satisfying_point_y)
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\beta$")
    plt.show()

def phase_boundary_degeneracy(coefs):
    satisfying_point_x=[]
    satisfying_point_y=[]
    Nv=12
    for alpha,beta in coefs:
        W=create_W_matrix3(alpha,beta)
        transferloop=diagonal_symmetry(W,Nv)
        m1s=np.sort(np.abs(np.linalg.eigvals(transferloop)))[::-1]
        m_temp=1-m1s[1]/m1s[0]
        if np.abs(m_temp)<=2E-6:
            satisfying_point_x.append(alpha)
            satisfying_point_y.append(beta)
    plt.scatter(satisfying_point_x,satisfying_point_y)
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\beta$")
    plt.show()

def test_phase_transition2Dplane(coefs):
    m1s=[]
    m2s=[]
    m3s=[]
    m4s=[]
    Nv=12
    fp1=open("./data/diagonal_symmetry{}.txt".format(Nv),"w")
    fp2=open("./data/diagonal_asymmetry{}.txt".format(Nv),"w")
    fp3=open("./data/offdiagonal_symmetry{}.txt".format(Nv),"w")
    fp4=open("./data/offdiagonal_asymmetry{}.txt".format(Nv),"w")
    for coef in coefs:
        alpha,beta=path1(coef)
        W=create_W_matrix3(alpha,beta)
        m1=np.sort(np.abs(np.linalg.eigvals(diagonal_symmetry(W,Nv))))[-1]
        m2=np.sort(np.abs(np.linalg.eigvals(diagonal_asymmetry(W,Nv))))[-1]
        m3=np.sort(np.abs(np.linalg.eigvals(offdiagonal_symmetry(W,Nv))))[-1]
        m4=np.sort(np.abs(np.linalg.eigvals(offdiagonal_asymmetry(W,Nv))))[-1]
        m1s.append(m1)
        m2s.append(m2)
        m3s.append(m3)
        m4s.append(m4)
        fp1.write("{}\t{}\n".format(coef,m1))
        fp2.write("{}\t{}\n".format(coef,m2))
        fp3.write("{}\t{}\n".format(coef,m3))
        fp4.write("{}\t{}\n".format(coef,m4))
    fp1.close()
    fp2.close()
    fp3.close()
    fp4.close()
    for i in range(len(m1s)):
        m4s[i]=m4s[i]/m1s[i]
        m3s[i]=m3s[i]/m1s[i]
        m2s[i]=m2s[i]/m1s[i]
        m1s[i]=m1s[i]/m1s[i]
    plt.plot(coefs,m1s,label=r"$|\gamma^{e\phi}_{e\phi}|$",marker='D')
    plt.plot(coefs,m2s,label=r"$|\gamma^{o\phi}_{o\phi}|$",marker='o')
    plt.plot(coefs,m3s,label=r"$|\gamma^{e\phi}_{e\phi+\pi}|$",marker='s')
    plt.plot(coefs,m4s,label=r"$|\gamma^{o\phi}_{o\phi+\pi}|$",marker='*')
    plt.legend()
    plt.show()


def test_alphabeta_gap():
    Nv=12
    coefs=np.linspace(0,1,20)
    theos=[]
    nums=[]
    for coef in coefs:
        W=create_W_matrix3(coef,coef)
        transferloop=diagonal_symmetry(W,Nv)
        m1s=np.sort(np.abs(np.linalg.eigvals(transferloop)))[::-1]
        m_temp=m1s[0]-m1s[1]
        nums.append(m_temp)
        #theos.append(4*(1-coef)*coef*(1-(2*coef-1)**(Nv-2)))
        if coef<=1/2:
            theos.append(2)
        theos.append(1+(2*coef-1)**Nv-np.abs(2*coef-1)*(1+(2*coef-1)**(Nv-2)))
    plt.plot(coefs,nums,label="Numerical",marker='D')
    plt.plot(coefs,theos,label="Theoritical",marker='o')
    plt.legend()
    plt.show()

def test_evenodd():
    Nv=11
    eigenplots1=[]
    eigenplots2=[]
    coefs=np.linspace(0,1,20)
    for coef in coefs:
        print(coef)
        alpha,beta=coef,coef
        W=create_W_matrix3(alpha,beta)
        loop_transfer=diagonal_phi_loop(W,Nv)
        loop_transfer2=diagonal_phi_loop(W,Nv-1)
        eigs=np.sort(np.abs(np.linalg.eigvals(loop_transfer)))[::-1][:32]
        print("The odds are:")
        print(eigs)
        eigs2=np.sort(np.abs(np.linalg.eigvals(loop_transfer2)))[::-1][:32]
        print("The even are:")
        print(eigs2)
#        print_eig=np.concatenate([eigs,eigs2])
        eigenplots1.append(eigs[0])
        eigenplots2.append(eigs2[0])
    plt.rcParams.update({'font.size': 25})
    plt.xlabel(r"$\alpha$")
    plt.plot(coefs,eigenplots1,label="Nv=11")
    plt.plot(coefs,eigenplots2,label="Nv=10")
    plt.legend()
    plt.show()

def test_diagonal_line():
    coefs=np.linspace(0,1,20)
    num1=[]
    num2=[]
    they1=[]
    they2=[]
    for coef in coefs:
        print(coef)
        W=create_W_matrix3(coef,coef)
        loop1=diagonal_phi_loop(W,11)
        loop2=diagonal_phi_loop(W,10)
        m1s=np.sort(np.abs(np.linalg.eigvals(loop1)))[::-1][:3]
        print(m1s)
        m1=m1s[2]
        m2=np.sort(np.abs(np.linalg.eigvals(loop2)))[-3]
        print(m1,m2)
        num1.append(m1)
        num2.append(m2)
        they1.append(np.sort(np.abs(np.array([(2*coef-1)**k+(2*coef-1)**(11-k) for k in range(12)])))[-3])
        they2.append(np.sort(np.abs(np.array([(2*coef-1)**k+(2*coef-1)**(10-k) for k in range(11)])))[-3])
    plt.plot(coefs,num1,label="Numerical Nv=11",marker='D')
    plt.plot(coefs,they1,label="Theoritical Nv=11",marker='o')
    plt.legend()
    plt.show()
    plt.plot(coefs,num2,label="Numerical Nv=10",marker='D')
    plt.plot(coefs,they2,label="Theoritical Nv=10",marker='o')
    plt.legend()
    plt.show()



def degeneracy_of_transfer2Dplane(coefs):
    Nv=11
    eigenplots1=[]
    eigenplots2=[]
    for coef in coefs:
        print(coef)
        alpha,beta=path1(coef)
        W=create_W_matrix3(alpha,beta)
        loop_transfer=diagonal_phi_loop(W,Nv)
        loop_transfer2=offdiagonal_phi_loop(W,Nv)
        eigs=np.sort(np.abs(np.linalg.eigvals(loop_transfer)))[::-1][:32]
        print("The 00/pipis are:")
        print(eigs)
        eigs2=np.sort(np.abs(np.linalg.eigvals(loop_transfer2)))[::-1][:32]
        print("The 0pi/pi0 are:")
        print(eigs2)
        if alpha==1 or beta==1:
            print("The diagonal error is {}".format(np.abs(1+coef**Nv-eigs[0])))
            print("The offdiagonal error is {}".format(np.abs(1-coef**Nv-eigs2[0])))
        if eigs[0]!=0:
            eigs2=eigs2/eigs[0]
        if eigs[0]!=0:
            eigs=eigs/eigs[0]
#        print_eig=np.concatenate([eigs,eigs2])
        eigenplots1.append(eigs)
        eigenplots2.append(eigs2)
    eigenplots1=np.array(eigenplots1)
    eigenplots2=np.array(eigenplots2)
    for i in range(32):
        plt.plot(coefs,eigenplots1[:,i])
    plt.show()
    for i in range(32):
        plt.plot(coefs,eigenplots2[:,i])
    plt.show()

def test_phase_transition(coefs):
    m1s=[]
    m2s=[]
    m3s=[]
    m4s=[]
    Nv=12
    fp1=open("./data/diagonal_symmetry{}.txt".format(Nv),"w")
    fp2=open("./data/diagonal_asymmetry{}.txt".format(Nv),"w")
    fp3=open("./data/offdiagonal_symmetry{}.txt".format(Nv),"w")
    fp4=open("./data/offdiagonal_asymmetry{}.txt".format(Nv),"w")
    for coef in coefs:
        W=create_W_matrix2(coef)
        m1=np.sort(np.abs(np.linalg.eigvals(diagonal_symmetry(W,Nv))))[-1]
        m2=np.sort(np.abs(np.linalg.eigvals(diagonal_asymmetry(W,Nv))))[-1]
        m3=np.sort(np.abs(np.linalg.eigvals(offdiagonal_symmetry(W,Nv))))[-1]
        m4=np.sort(np.abs(np.linalg.eigvals(offdiagonal_asymmetry(W,Nv))))[-1]
        m1s.append(m1)
        m2s.append(m2)
        m3s.append(m3)
        m4s.append(m4)
        fp1.write("{}\t{}\n".format(coef,m1))
        fp2.write("{}\t{}\n".format(coef,m2))
        fp3.write("{}\t{}\n".format(coef,m3))
        fp4.write("{}\t{}\n".format(coef,m4))
    fp1.close()
    fp2.close()
    fp3.close()
    fp4.close()
    for i in range(len(m1s)):
        m4s[i]=m4s[i]/m1s[i]
        m3s[i]=m3s[i]/m1s[i]
        m2s[i]=m2s[i]/m1s[i]
        m1s[i]=m1s[i]/m1s[i]
    plt.rcParams.update({'font.size': 25})
    plt.plot(coefs,m1s,label=r"$|\gamma^{e\phi}_{e\phi}|$",marker='D')
    plt.plot(coefs,m2s,label=r"$|\gamma^{o\phi}_{o\phi}|$",marker='o')
    plt.plot(coefs,m3s,label=r"$|\gamma^{e\phi}_{e\phi+\pi}|$",marker='s')
    plt.plot(coefs,m4s,label=r"$|\gamma^{o\phi}_{o\phi+\pi}|$",marker='*')
    plt.xlabel(r"$\theta$")
    plt.legend()
    plt.show()

def degeneracy_of_transfer(coefs):
    Nv=12
    eigenplots=[]
    for coef in coefs:
        W=create_W_matrix2(coef)
        loop_transfer=diagonal_phi_loop(W,Nv)
        loop_transfer2=offdiagonal_phi_loop(W,Nv)
        eigs=np.sort(np.abs(np.linalg.eigvals(loop_transfer)))[::-1][:32]
#        eigs=np.sort(np.linalg.eigvals(loop_transfer))[:32]
        print("The 00/pipis are:")
        print(eigs)
        if eigs[0]!=0:
            eigs=eigs/eigs[0]
#        eigs2=np.sort(np.abs(np.linalg.eigvals(loop_transfer2)))[::-1][:32]
#        if eigs[0]!=0:
#            eigs2=eigs2/eigs[0]
#        print("The 0pi/pi0 are:")
#        print(eigs2)
#        print_eig=np.concatenate([eigs,eigs2])
        print_eig=eigs
        eigenplots.append(print_eig)
    eigenplots=np.array(eigenplots)
    plt.rcParams.update({'font.size': 25})
    plt.xlabel(r"$\theta$")
    for i in range(32):
        plt.plot(coefs,eigenplots[:,i])
    plt.show()


def Gapclosesrate_differentblocks(theta=0.05):
    fp=open("./data/gapcloserateDifferentBlock{}.txt".format(theta),"w")
    Nvs=[]
    ms=[]
    W=create_W_matrix2(theta)
    for Nv in range(2,13):
        m1=np.sort(np.abs(np.linalg.eigvals(diagonal_symmetry(W,Nv))))[-1]
        m2=np.sort(np.abs(np.linalg.eigvals(offdiagonal_symmetry(W,Nv))))[-1]
        print(m1,m2)
        m_temp=1-m2/m1
        fp.write("{}\t{}\n".format(Nv,m_temp))
        Nvs.append(Nv)
        ms.append(m_temp)
    fp.close()
    plt.plot(Nvs,ms)
    plt.show()

def Gapcloserate_sameblocks(theta=0.05):
    fp=open("./data/gapcloserateSameBlock{}.txt".format(theta),"w")
    Nvs=[]
    ms=[]
    W=create_W_matrix2(theta)
    for Nv in range(2,13):
        transferloop=diagonal_symmetry(W,Nv)
        m1s=np.sort(np.abs(np.linalg.eigvals(transferloop)))[::-1]
        print(m1s[:50])
        m_temp=1-m1s[1]/m1s[0]
        fp.write("{}\t{}\n".format(Nv,m_temp))
        Nvs.append(Nv)
        ms.append(m_temp)
    fp.close()
    plt.plot(Nvs,ms)
    plt.show()


def sameblock_analytic_verification(coefs):
    Nv=12
    plotdatas=[]
    theoritical_datas=[]
    for coef in coefs:
        W=create_W_matrix2(coef)
        transferloop=diagonal_symmetry(W,Nv)
        m1s=np.sort(np.abs(np.linalg.eigvals(transferloop)))[::-1]
        m_temp=1-m1s[1]/m1s[0]
        plotdatas.append(m_temp)
        theoritical_datas.append(np.sin(2*coef)**2)
    plt.plot(coefs,plotdatas,label="Numerical",marker='D')
    plt.plot(coefs,theoritical_datas,label="Theoritical",marker='o')
    plt.legend()
    plt.show()

def differentblock_analytic_verification(coef):
    Nv=12
    plotdatas=[]
    theoritical_datas=[]
    for coef in coefs:
        W=create_W_matrix2(coef)
        m1=np.sort(np.abs(np.linalg.eigvals(diagonal_symmetry(W,Nv))))[-1]
        m2=np.sort(np.abs(np.linalg.eigvals(offdiagonal_symmetry(W,Nv))))[-1]
        m_temp=1-m2/m1
        plotdatas.append(m_temp)
        theoritical_datas.append(1-np.cos(2*coef))
    plt.plot(coefs,plotdatas,label="Numerical",marker='D')
    plt.plot(coefs,theoritical_datas,label="Theoritical",marker='o')
    plt.legend()
    plt.show()

def perturbed_around_alpha1():
    numes=[]
    theos=[]
    Nv=12
    beta=0.0
    alphas=[]
    for alpha in np.linspace(1,0.95,10):
        W=create_W_matrix3(alpha,beta)
        m1=np.sort(np.abs(np.linalg.eigvals(diagonal_symmetry(W,Nv))))[-1]
        m2=np.sort(np.abs(np.linalg.eigvals(offdiagonal_symmetry(W,Nv))))[-1]
        print("alpha={},beta={}".format(alpha,beta))
        print("Numerical result is {}\t{}",m1,m2)
        m_temp=1-m2/m1
        numes.append(m_temp)
        delta=1-alpha
        identi=beta**Nv+(1-delta)**Nv
        first_term=identi
        for j in range(1,Nv):
            first_term += Nv*delta*beta**(j-1)*(1-beta)*(1-delta)**(Nv-j-1)
        second_term=-beta**Nv+(1-delta)**Nv
        for i in range(1,Nv+1):
            for j in range(i+1,Nv+1):
                second_term+=alpha**(Nv-2)*delta*(1-beta)*((beta/alpha)**(j-i-1)-(beta/alpha)**(Nv-j+i-1))
        print("Theoritical result is {}\t{}".format(first_term,second_term))
        theos.append(1-second_term/first_term)
        alphas.append(alpha)
    plt.plot(alphas,numes,label="Numerical",marker='D')
    plt.plot(alphas,theos,label="Theoritical",marker='o')
    plt.legend()
    plt.show()


if __name__=='__main__':
    #perturbed_around_alpha1()
    coefs=np.linspace(0,np.pi/4,10)
    #coefs=np.array([0.05,0.2,0.4,0.5,0.7,0.72,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.88,0.9,0.95,1.0])
    #test_phase_transition(coefs)
    degeneracy_of_transfer(coefs)
    #sameblock_analytic_verification(coefs)
    #differentblock_analytic_verification(coefs)
    #test_phase_transition2Dplane(coefs)
    #degeneracy_of_transfer2Dplane(coefs)
    coefs=[]
    for alpha in np.linspace(0,1.0,10):
        for beta in np.linspace(0,1.0,10):
            coefs.append((alpha,beta))
    phase_boundary_overlap(coefs)
    #coefs=[]
    #for alpha in np.linspace(0,1,10):
    #    for beta in np.linspace(0,1,10):
    #        coefs.append((alpha,beta))
    #phase_boundary_overlap(coefs)
