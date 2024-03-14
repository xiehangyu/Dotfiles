import matplotlib.pyplot as plt
import numpy as np
def energy_entropy_plot():
    onsite_random=input("onsite_random:\n")
    L=input("L:\n")
    alpha=input("alpha:\n")
    modify_name=input("modify_name:\n")
    D=input("D:\n")
    fp=open("./plot_data/{}EnergyD{}L{}alpha{}onsiterandom{}.txt".format(modify_name,D,L,alpha,onsite_random),'r')
    data=fp.read().split()
    data=[float(i) for i in data]
    y=np.zeros(len(data))
    plt.plot(y,data,'o')
    plt.show()
    fp=open("./plot_data/{}EntropyD{}L{}alpha{}onsiterandom{}.txt".format(modify_name,D,L,alpha,onsite_random),'r')
    data=fp.read().split()
    data=[float(i) for i in data]
    x=data[0::2]
    y=data[1::2]
    plt.plot(x,y)
    plt.show()


if __name__=='__main__':
    energy_entropy_plot()