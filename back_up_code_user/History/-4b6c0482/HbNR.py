import matplotlib.pyplot as plt
import numpy as np

def temp1():
    plt.rcParams.update({'font.size': 30})
    h=0.0
    t=0.02
    tend=60
    randomin=0
    Korder=3
    for alpha in [0.4]:
        plt.figure(figsize=(12,8))
        if alpha==0.0:
            GSE=24
        else:
            GSE=22
        for N in [20,40,60]:
            if N==20 or N==40:
                tend=50
                GSE=160
                Korder=5
            else:
                tend=60
                GSE=22
                Korder=3
            filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            if alpha==0.0 and N==60:
                GSE_temp=160
                tend_temp=50
                Korder_temp=5
                filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend_temp,GSE_temp,Korder_temp,randomin,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.plot(ts/4,Entropys,label='N={}'.format(N))
        plt.xlabel(r't')
        plt.xlim([0,12.5])
        plt.ylabel(r'$S_A$')
        plt.legend()
        if alpha==-100:
            plt.title('alpha=inf,h={}'.format(h*2))
        else:
            plt.title('alpha={},h={}'.format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/alpha{}h{}.png".format(alpha,h*2))
        plt.show()

def temp2():
    N=60
    t=0.02
    tend=60
    GSE=22
    Korder=3
    randomin=0
    h=0.1
    plt.rcParams.update({'font.size': 30})
    plt.figure(figsize=(12,8))
    for alpha in [0.0,0.1,0.2,0.4,0.8,1.0,-100]:
        filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend,GSE,Korder,randomin,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        if alpha==-100:
            plt.plot(ts/4,Entropys,label='alpha=inf')
        else:
            plt.plot(ts/4,Entropys,label='alpha={}'.format(alpha))
    plt.xlim([0,10])
    plt.xlabel('t')
    plt.ylabel(r'$S_A$')
    plt.legend()
    plt.title('h={},N={}'.format(h*2,N))
    #plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/h{}N{}.png".format(h*2,N))
    plt.show()

def temp3():
    h=0.0
    N=60
    t=0.02
    tend=60
    GSE=22
    Korder=3
    randomin=0
    plt.rcParams.update({'font.size': 30})
    plt.figure(figsize=(12,8))
    for alpha in [0.0,0.1,0.2,0.4,0.8,1.0,-100]:
    #for alpha in [0.0,0.1,0.2,0.4,0.8,1.0,2.0,-100]:
        filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
        if alpha==0.0 or alpha==0.4 or alpha==2.0:
            GSE_temp=160
            tend_temp=50
            Korder_temp=5
            filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend_temp,GSE_temp,Korder_temp,randomin,h)
        print(filename)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        if alpha==-100:
        #    plt.plot(ts/4,np.exp(Entropys),label='alpha=inf')
            plt.semilogx(ts/4,Entropys,label='alpha=inf')
        else:
           # plt.plot(ts/4,np.exp(Entropys),label='alpha={}'.format(alpha))
            plt.semilogx(ts/4,Entropys,label='alpha={}'.format(alpha))
    plt.xlabel('t')
    plt.ylabel(r'$\exp{S_A}$')
    plt.legend()
    plt.title('h={},N={}'.format(h*2,N))
    #plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/Exph{}N{}.png".format(h*2,N))
    plt.show()

def temp4():
    GSE=25
    Korder=3
    randomin=0
    t=0.2
    tend=60
    h=0.5
    plt.rcParams.update({'font.size': 30})
    for alpha in [0.1,0.7]:
        plt.figure(figsize=(12,8))
        for N in [20,40,60,80]:
            if h==0.2:
                tend=60
                GSE=25
                t=0.02
            else:
                t=0.02
                tend=60
                GSE=25
            filename='NewXXOBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            print(filename)
            #filename="XXexactN_{}t_{}tend_{}h_{}.txt".format(N,t,tend,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            #plt.semilogx(ts,Entropys,label='N={}'.format(N))
            plt.plot(ts,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlim([0,5])
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XXPeakcompared_a_alpha{}h{}.png".format(alpha,h))
        plt.show()


def temp5():
    GSE=25
    Korder=3
    randomin=0
    t=0.2
    tend=60
    N=60
    alpha=-100
    plt.rcParams.update({'font.size': 30})
    for h in [0.01,0.02,0.04,0.08,0.12,0.16]:
        filename='NewXXOBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
        print(filename)
        #filename="XXexactN_{}t_{}tend_{}h_{}.txt".format(N,t,tend,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        #plt.semilogx(ts,Entropys,label='N={}'.format(N))
        plt.plot(ts,Entropys,label='h={}'.format(h))
    plt.rcParams.update({'font.size': 30})
    plt.xlim([0,5])
    plt.xlabel('t')
    plt.ylabel(r"$S_A$")
    plt.legend()
    if alpha==-100:
        plt.title(r"$\alpha={},N={}$".format("inf",N))
    else:
        plt.title(r"$\alpha={},N={}$".format(alpha,N))
    if alpha==-100:
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XXPeakcompared_a_alphainfN{}.png".format(N))
    else:
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XXPeakcompared_a_alpha{}N{}.png".format(alpha,N))
    plt.show()

if __name__ == '__main__':
    temp5()