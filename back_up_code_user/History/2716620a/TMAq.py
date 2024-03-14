import numpy as np
import matplotlib.pyplot as plt

def draw_distance(N,alpha):
    fp=open("distancelongrange_N_"+str(N)+"_alpha_"+str(alpha)+".txt","r")
    data=fp.read().split()
    times=[float(i) for i in data[0::4]]
    sqrttracedistance=[complex(i) for i in data[1::4]]
    sqrthilbertdistance=[complex(i) for i in data[2::4]]
    tracedistance=[complex(i) for i in data[3::4]]
    plt.plot(times,sqrttracedistance,label="sqrttracedistance")
    plt.plot(times,sqrthilbertdistance,label="sqrthilbertdistance")
    plt.plot(times,tracedistance,label="tracedistance")
    plt.xlabel("time")
    plt.ylabel("distance")
    plt.legend()
    plt.show()
    fp=open("distancenearest_N_"+str(N)+".txt","r")
    data=fp.read().split()
    times=[float(i) for i in data[0::4]]
    sqrttracedistance=[complex(i) for i in data[1::4]]
    sqrthilbertdistance=[complex(i) for i in data[2::4]]
    tracedistance=[complex(i) for i in data[3::4]]
    plt.plot(times,sqrttracedistance,label="sqrttracedistance")
    plt.plot(times,sqrthilbertdistance,label="sqrthilbertdistance")
    plt.plot(times,tracedistance,label="tracedistance")
    plt.xlabel("time")
    plt.ylabel("distance")
    plt.legend()
    plt.show()

def draw_ratio(N,alpha):
    fp=open("distancelongrange_N_"+str(N)+"_alpha_"+str(alpha)+".txt","r")
    data=fp.read().split()
    times=[float(i) for i in data[0::4]]
    sqrttracedistance=[complex(i) for i in data[1::4]]
    sqrthilbertdistance=[complex(i) for i in data[2::4]]
    tracedistance=[complex(i) for i in data[3::4]]
    ratios=[]
    for i in range(len(sqrttracedistance)):
        if sqrthilbertdistance[i]==0:
            ratios.append(0)
        else:
            ratios.append(sqrttracedistance[i]/sqrthilbertdistance[i])
    plt.plot(times,ratios,label="ratio")
    plt.xlabel("time")
    plt.ylabel("ratio")
    plt.legend()
    plt.show()
    fp=open("distancenearest_N_"+str(N)+".txt","r")
    data=fp.read().split()
    times=[float(i) for i in data[0::4]]
    sqrttracedistance=[complex(i) for i in data[1::4]]
    sqrthilbertdistance=[complex(i) for i in data[2::4]]
    tracedistance=[complex(i) for i in data[3::4]]
    ratios=[]
    for i in range(len(sqrttracedistance)):
        if sqrthilbertdistance[i]==0:
            ratios.append(0)
        else:
            ratios.append(sqrttracedistance[i]/sqrthilbertdistance[i])
    plt.plot(times,ratios,label="ratio")
    plt.xlabel("time")
    plt.ylabel("ratio")
    plt.legend()
    plt.show()

def draw_fitting():
    NAlist=[]
    ratioslongrange=[]
    ratiosnearest=[]
    for i in range(10,24,2):
        NAlist.append(i/2)
        fp=open("distancelongrange_N_"+str(i)+"_alpha_0.1.txt","r")
        data=fp.read().split()
        sqrttracedistance=[complex(i) for i in data[1::4]]
        sqrthilbertdistance=[complex(i) for i in data[2::4]]
        ratioslongrange.append(sqrttracedistance[-1]/sqrthilbertdistance[-1])
        fp=open("distancenearest_N_"+str(i)+".txt","r")
        sqrttracedistance=[complex(i) for i in data[1::4]]
        sqrthilbertdistance=[complex(i) for i in data[2::4]]
        ratiosnearest.append(sqrttracedistance[-1]/sqrthilbertdistance[-1])
    plt.plot(NAlist,ratioslongrange,label="longrange")
    plt.xlabel(r"$N_A$")
    plt.ylabel(r"$\frac{d_{tr}}{d_{H}}$")
    plt.show()
    plt.plot(NAlist,ratiosnearest,label="nearest")
    plt.xlabel(r"$N_A$")
    plt.ylabel(r"$\frac{d_{tr}}{d_{H}}$")
    plt.show()


if __name__=='__main__':
    decision=input("1 for distance and drawing\n2 for ratio drawing\n")
    N=input("input N\n")
    alpha=input("input alpha\n")
    if decision=='1':
        draw_distance(N,alpha)
    if decision=='2':
        draw_ratio(N,alpha)
    if decision=='3':
        draw_fitting()