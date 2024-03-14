import matplotlib.pyplot as plt
import numpy as np
def energy_entropy_plot():
    onsite_random=input("onsite_random")
    L=input("L")
    alpha=input("alpha")
    modify_name=input("modify_name")
    D=input("D")
    fp=open("{}EnergyD{}L{}alpha{}onsiterandom{}.txt".format(modify_name,D,L,alpha,onsite_random),'r')
    data=fp.read().split()
    data=[float(i) for i in data]
    y=np.zeros(len(data))
    plt.plot(data,y,'o')
    plt.show()
    fp=open("{}EntropyD{}L{}alpha{}onsiterandom{}.txt".format(modify_name,D,L,alpha,onsite_random),'r')
    data=fp.read().split()
    data=[float(i) for i in data]
    x=data(0::2)
    y=data(1::2)
    plt.plot(x,y)
    plt.show()


if __name__=='__main__':
    energy_entropy_plot()