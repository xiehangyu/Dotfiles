import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap
def drawing_entanglement(N,alpha,t,tend,GSE,Korder,randomin,J):
    filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,J)
    fp=open(filename,'r')
    data=fp.read().split()
    ts=np.array(data[0::3],dtype=float)
    Entropys=np.array(data[2::3],dtype=float)
    plt.plot(ts,Entropys)
    plt.show()



def drawing_correlators(N,alpha,t,tend,GSE,Korder,J):
    filename='QuenchedCorrelatorsN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,J)
    fp=open(filename,'r')
    data=fp.read().split()
    print("The energy is {}".format(data[0]))
    ts=np.array(data[1::7],dtype=float)
    XXs=np.array(data[2::7],dtype=float)
    for i in range(np.size(ts)):
        XXs[i]*=(-1)**(i+1)
    plt.plot(ts,XXs)
    plt.show()

def uniform_in_one_picture():
    cmap = get_cmap('viridis')
    for h in range(11):
        color = cmap(h / (11 - 1)) 
        filename='EntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(40,0.1,0.2,120,40,4,0,h*0.1)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        plt.plot(ts,Entropys,color=color)
    plt.legend([str(h*0.1) for h in range(11)])
    plt.show()


def uniform_in_one_picture1():
    alpha=0
    GSE=24
    Korder=3
    N=60
    randomin=0
    t=0.02
    tend=60
    plt.figure(figsize=(12,8))
    for h in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]:
        if alpha==0:
            if h==0.0 or h==0.4 or h==0.8 or h==0.2:
                GSE=160
                Korder=5
                tend=50
        filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        plt.plot(ts/4,Entropys,label='h={}'.format(h*2))
        GSE=24
        Korder=3
        tend=60
    plt.xlabel('t')
    plt.ylabel(r'$S_A$')
    plt.legend()
    plt.title('alpha={},N={}'.format(alpha,N))
    plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/alpha{}N{}.png".format(alpha,N))
    plt.show()
    for alpha in [0.1,1,-100]:
        plt.figure(figsize=(12,8))
        GSE=22
        for h in [0.0,0.2,0.4,0.6,0.8]:
            filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend,GSE,Korder,randomin,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.plot(ts/4,Entropys,label='h={}'.format(h*2))
        plt.xlabel('t')
        plt.ylabel(r'$S_A$')
        plt.legend()
        if alpha==-100:
            plt.title('alpha=inf,N={}'.format(N))
        else:
            plt.title('alpha={},N={}'.format(alpha,N))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/alpha{}N{}.png".format(alpha,N))
        plt.show()

def compared_with_PRRxx():
    alpha=0
    GSE=22
    Korder=3
    randomin=0
    t=0.2
    tend=120
    for h in [0.2,0.4,0.5,0.8]:
        plt.figure(figsize=(12,8))
        for N in [50,100,200,400,800]:
            if h>0.3 and N==800:
                continue
            #filename='XXOBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            filename="XXexactN_{}t_{}tend_{}h_{}.txt".format(N,t,tend,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::2],dtype=float)
            Entropys=np.array(data[1::2],dtype=float)
            #plt.semilogx(ts,Entropys,label='N={}'.format(N))
            plt.plot(ts,Entropys,label='N={}'.format(N))
        if h==0.2 or h==0.4:
            filename="XXspinwave_N{}_alpha{}_t{}_tend{}_h{}.txt".format(100,0.0,t,60.0,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::2],dtype=float)
            Entropys=np.array(data[1::2],dtype=float)
            plt.semilogx(ts,Entropys,label="SW")
        if h==0.4:
            plt.xlim([0,15])
        else:
            plt.xlim([0,60])
        plt.rcParams.update({'font.size': 30})
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XXexactcompared_a_alpha{}h{}.png".format(alpha,h))
        plt.show()
    t=0.02
    tend=120
    GSE=25
    for h in [0.2,0.5]:
        if h==0.5:
            t=0.02
            tend=120
            GSE=25
            break

        plt.figure(figsize=(12,8))
        for N in [50,100,200,400,800]:
            if h==0.2:
                tend=60
                GSE=25
                t=0.02
            else:
                t=0.02
                tend=120
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
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XXcompared_a_alpha{}h{}.png".format(alpha,h))
        plt.show()
    for h in [0.2,0.5]:
        if h==0.5:
            t=0.02
            tend=120
            GSE=25
        plt.figure(figsize=(12,8))
        for N in [50,100,200,400,800]:
            if h==0.2:
                tend=60
                GSE=25
                t=0.2
            else:
                t=0.02
                tend=120
                GSE=25
            filename='XXOBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            print(filename)
            #filename="XXexactN_{}t_{}tend_{}h_{}.txt".format(N,t,tend,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.semilogx(ts/4,Entropys,label='N={}'.format(N))
            #plt.plot(ts,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XXcompared_a_alpha{}h{}.png".format(alpha,h*2))
        plt.show()
    for h in [0.2,0.5]:
        if h==0.5:
            t=0.02
            tend=120
            GSE=25
            break
        plt.figure(figsize=(12,8))
        for N in [50,100,200,400,800]:
            if h==0.2:
                tend=60
                GSE=0
                t=0.02
            else:
                t=0.02
                tend=120
                GSE=25
                break
            filename='XXOBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            print(filename)
            #filename="XXexactN_{}t_{}tend_{}h_{}.txt".format(N,t,tend,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.semilogx(ts/4,Entropys,label='N={}'.format(N))
            #plt.plot(ts,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XX2compared_a_alpha{}h{}.png".format(alpha,h*2))
        plt.show()
    for h in [0.2,0.5]:
        if h==0.5:
            t=0.02
            tend=120
            GSE=25
            break
        plt.figure(figsize=(12,8))
        for N in [50,100,200,400,800]:
            if h==0.2:
                tend=60
                GSE=250
                t=0.02
            else:
                t=0.02
                tend=120
                GSE=25
                break
            filename='XXOBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            print(filename)
            #filename="XXexactN_{}t_{}tend_{}h_{}.txt".format(N,t,tend,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.semilogx(ts/4,Entropys,label='N={}'.format(N))
            #plt.plot(ts,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XX3compared_a_alpha{}h{}.png".format(alpha,h*2))
        plt.show()
    GSE=22
    alpha=0.1
    tend=60
    for h in [0.0,0.2,0.5,2.0]:
        plt.figure(figsize=(12,8))
        for N in [20,40,60,80]:
            filename='XXOBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.semilogx(ts/4,Entropys,label='N={}'.format(N))
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.rcParams.update({'font.size': 30})
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XXcompared_c_alpha{}h{}.png".format(alpha,h*2))
        plt.show()
    alpha=0.7
    GSE=22
    for h in [0.0,0.2,0.5,2.0]:
        plt.figure(figsize=(12,8))
        for N in [20,40,60,80]:
            filename='XXOBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            print(filename)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.semilogx(ts/4,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/XXcompared_d_alpha{}h{}.png".format(alpha,h*2))
        plt.show()
def compared_with_PRR():
    alpha=0
    t=0.2
    tend=120
    plt.rcParams.update({'font.size': 30})
    for h in [0.2,0.4,0.5,0.8]:
        plt.figure(figsize=(12,8))
        for N in [50,100,200,400]:
            if h>0.3 and N==800:
                continue
            #filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            filename="exactN_{}t_{}tend_{}h_{}.txt".format(N,t,tend,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::2],dtype=float)
            Entropys=np.array(data[1::2],dtype=float)
            plt.semilogx(ts,Entropys,label='N={}'.format(N))
            #plt.plot(ts,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlim([0,60])
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/exactcompared_a_semilog_alpha{}h{}.png".format(alpha,h))
        plt.show()
    t=0.02
    tend=60
    Korder=3
    GSE=25
    randomin=0.0
    for h in [0.2]:
        plt.figure(figsize=(12,8))
        for N in [50,100,200,400,800]:
            filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            print(filename)
            #filename="exactN_{}t_{}tend_{}h_{}.txt".format(N,t,tend,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.semilogx(ts,Entropys,label='N={}'.format(N))
            #plt.plot(ts,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r"$\alpha={},h={}$".format(alpha,h))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/compared_a_alpha{}h{}.png".format(alpha,h))
        plt.show()
    GSE=24
    Korder=3
    randomin=0
    t=0.02
    tend=60
    h=0.1
    plt.rcParams.update({'font.size': 30})
    plt.figure(figsize=(12,8))
    for N in [20,30,40,60,80]:
        filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        #plt.semilogx(ts,Entropys,label='N={}'.format(N))
        plt.plot(ts/4,Entropys,label='N={}'.format(N))
    plt.xlabel('t')
    plt.ylabel(r"$S_A$")
    plt.legend()
    plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
    plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/compared_a_alpha{}h{}.png".format(alpha,h*2))
    plt.show()
    h=0.2
    plt.figure(figsize=(12,8))
    for N in [20,30,40,60,80]:
        if N!=60:
            filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
        else:
            filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,50,160,5,randomin,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        plt.plot(ts/4,Entropys,label='N={}'.format(N))
    plt.rcParams.update({'font.size': 30})
    plt.xlabel('t')
    plt.ylabel(r"$S_A$")
    plt.legend()
    plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
    plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/compared_b_alpha{}h{}.png".format(alpha,h*2))
    plt.show()
    h=0.5
    GSE=22
    plt.figure(figsize=(12,8))
    for N in [20,40,60,80]:
        filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        plt.plot(ts/4,Entropys,label='N={}'.format(N))
    plt.rcParams.update({'font.size': 30})
    plt.xlabel('t')
    plt.ylabel(r"$S_A$")
    plt.legend()
    plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
    plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/compared_b_alpha{}h{}.png".format(alpha,h*2))
    plt.show()
    h=0.8
    plt.figure(figsize=(12,8))
    for N in [20,40,60,80]:
        filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        plt.plot(ts/4,Entropys,label='N={}'.format(N))
    plt.rcParams.update({'font.size': 30})
    plt.xlabel('t')
    plt.ylabel(r"$S_A$")
    plt.legend()
    plt.title(r"$\alpha={},h={}$".format(alpha,h*2))
    plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/compared_b_alpha{}h{}.png".format(alpha,h*2))
    plt.show()
    alpha=0.1
    t=0.02
    tend=60
    h=0.2
    GSE=22
    for h in[0.2,0.8]:
        plt.figure(figsize=(12,8))
        for N in [20,40,60,80]:
            filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.semilogx(ts,Entropys,label='N={}'.format(N))
            #plt.plot(ts/4,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlabel('t')
        plt.ylabel(r'$S_A$')
        plt.legend()
        plt.title(r'$\alpha={},h={}$'.format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/compared_c_semilog_alpha{}h{}.png".format(alpha,h*2))
        plt.show()
    alpha=0.7
    t=0.02
    tend=60
    h=0.2
    GSE=22
    for h in [0.2,0.8]:
        plt.figure(figsize=(12,8))
        for N in [20,40,60,80]:
            filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
            print(filename)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.semilogx(ts/4,Entropys,label='N={}'.format(N))
        plt.rcParams.update({'font.size': 30})
        plt.xlabel('t')
        plt.ylabel(r"$S_A$")
        plt.legend()
        plt.title(r'$\alpha={},h={}$'.format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/compared_d_semilog_alpha{}h{}.png".format(alpha,h*2))
        plt.show()


def uniform_in_one_picture3():
    h=0.2
    t=0.02
    tend=60
    randomin=0
    Korder=3
    for alpha in [0.0,0.1,1.0,-100]:
        plt.figure(figsize=(12,8))
        if alpha==0.0:
            GSE=24
        else:
            GSE=22
        for N in [20,40,60,80]:

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
            plt.plot(ts/4,np.exp(Entropys),label='N={}'.format(N))
        plt.xlabel(r't')
        plt.ylabel(r'$\exp{S_A}$')
        plt.legend()
        if alpha==-100:
            plt.title('alpha=inf,h={}'.format(h*2))
        else:
            plt.title('alpha={},h={}'.format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/Exp_alpha{}h{}.png".format(alpha,h*2))
        plt.show()
    h=0.0
    GSE=22
    for alpha in [0.0,0.1,1.0,-100]:
        plt.figure(figsize=(12,8))
        for N in [20,40,60,80]:
            if alpha==0.0 and N!=80:
                GSE_temp=160
                tend_temp=50
                Korder_temp=5
                filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend_temp,GSE_temp,Korder_temp,randomin,h)
            else:
                filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend,GSE,Korder,randomin,h)
            fp=open(filename,'r')
            data=fp.read().split()
            ts=np.array(data[0::3],dtype=float)
            Entropys=np.array(data[2::3],dtype=float)
            plt.plot(ts/4,np.exp(Entropys),label='N={}'.format(N))
        plt.xlabel(r't')
        plt.ylabel(r'$\exp{S_A}$')
        plt.legend()
        if alpha==-100:
            plt.title('alpha=inf,h={}'.format(h*2))
        else:
            plt.title('alpha={},h={}'.format(alpha,h*2))
        plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/Exp_alpha{}h{}.png".format(alpha,h*2))
        plt.show()
    



def uniform_in_one_picture2():
    h=0.0
    N=60
    t=0.02
    tend=60
    GSE=22
    Korder=3
    randomin=0
    plt.figure(figsize=(12,8))
    for alpha in [0.0,0.1,0.4,0.8,1.0,2.0,-100]:
        filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,h)
        if alpha==0.0 or alpha==0.4 or alpha==0.8 or alpha==2.0:
            GSE_temp=160
            tend_temp=50
            Korder_temp=5
            filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend_temp,GSE_temp,Korder_temp,randomin,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        if alpha==-100:
            plt.plot(ts/4,Entropys,label='alpha=inf')
        else:
            plt.plot(ts/4,Entropys,label='alpha={}'.format(alpha))
    plt.xlabel('t')
    plt.ylabel(r'$S_A$')
    plt.legend()
    plt.title('h={},N={}'.format(h*2,N))
    plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/h{}N{}.png".format(h*2,N))
    plt.show()
    h=0.2
    plt.figure(figsize=(12,8))
    for alpha in [0.0,0.1,0.5,1.0,2.0,-100]:
        filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend,GSE,Korder,randomin,h)
        if alpha==0.0 or alpha==0.5 or alpha==2.0:
            GSE_temp=160
            tend_temp=50
            Korder_temp=5
            filename="OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt".format(N,alpha,t,tend_temp,GSE_temp,Korder_temp,randomin,h)
        fp=open(filename,'r')
        data=fp.read().split()
        ts=np.array(data[0::3],dtype=float)
        Entropys=np.array(data[2::3],dtype=float)
        if alpha==-100:
            plt.plot(ts/4,Entropys,label='alpha=inf')
        else:
            plt.plot(ts/4,Entropys,label='alpha={}'.format(alpha))
    plt.xlabel('t')
    plt.ylabel(r'$S_A$')
    plt.legend()
    plt.title('h={},N={}'.format(h*2,N))
    plt.savefig("/home/xiehangyu/yxh/study/mpq/Long_range_area_law/script/figs/h{}N{}.png".format(h*2,N))
    plt.show()
    


def drawing_eigenvalues(N,alpha,t,tend,GSE,Korder,randomin,J):
    filename='OBCEntropyresultN{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.txt'.format(N,alpha,t,tend,GSE,Korder,randomin,J)
    fp=open(filename,'r')
    data=fp.read().split()
    ts=np.array(data[0::3],dtype=float)
    Entropys=np.array(data[1::3],dtype=float)
    plt.plot(ts,Entropys)
    plt.show()

if __name__=='__main__':
    #compared_with_PRR()
    compared_with_PRRxx()
    uniform_in_one_picture2()
    uniform_in_one_picture1()
    uniform_in_one_picture3()
    for alpha in [0.0,1.0,2.0,-100]:
        for h in [0.0,0.8,1.6]:
            for N in [30,60]:
                print(N,alpha,h)
                drawing_entanglement(N,alpha,0.02,50,160,5,0,h)
    for alpha in [0.0,1.0,2.0,-100]:
        for h in [0.0,0.8,1.6]:
            for N in [30,60]:
                print(N,alpha,h)
                drawing_eigenvalues(N,alpha,0.02,50,160,5,0,h)
    #uniform_in_one_picture()
 #   for alpha in range(21):
 #       print(alpha/10)
 #       drawing_correlators(40,alpha*0.1,0.3,240,160,5,0)


    