import numpy as np
import matplotlib.pyplot as plt

def draw_distance(N,alpha):
    fp=open("distancelongrange_N_"+str(N)+"_alpha_"+str(alpha)+".txt","r")
    data=fp.read().split()
    times=[float(i) for i in data[0::4]]
    sqrttracedistance=[float(i) for i in data[1::4]]
    sqrthilbertdistance=[float(i) for i in data[2::4]]
    tracedistance=[float(i) for i in data[3::4]]
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
    sqrttracedistance=[float(i) for i in data[1::4]]
    sqrthilbertdistance=[float(i) for i in data[2::4]]
    tracedistance=[float(i) for i in data[3::4]]
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
    sqrttracedistance=[float(i) for i in data[1::4]]
    sqrthilbertdistance=[float(i) for i in data[2::4]]
    tracedistance=[float(i) for i in data[3::4]]
    ratios=[]
    for i in range(len(sqrttracedistance)):
        ratios.append(sqrttracedistance[i]/sqrthilbertdistance[i])
    plt.plot(times,ratios,label="ratio")
    plt.xlabel("time")
    plt.ylabel("ratio")
    plt.legend()
    plt.show()
    fp=open("distancenearest_N_"+str(N)+".txt","r")
    data=fp.read().split()
    times=[float(i) for i in data[0::4]]
    sqrttracedistance=[float(i) for i in data[1::4]]
    sqrthilbertdistance=[float(i) for i in data[2::4]]
    tracedistance=[float(i) for i in data[3::4]]
    ratios=[]
    for i in range(len(sqrttracedistance)):
        ratios.append(sqrttracedistance[i]/sqrthilbertdistance[i])
    plt.plot(times,ratios,label="ratio")
    plt.xlabel("time")
    plt.ylabel("ratio")
    plt.legend()
    plt.show()



if __name__=='__main__':
    decision=input("1 for distance and drawing\n 2 for ratio drawing\n")
    N=input("input N\n")
    alpha=input("input alpha\n")
    if decision==1:
        draw_distance(N,alpha)