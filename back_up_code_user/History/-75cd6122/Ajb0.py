import numpy as np
import matplotlib.pyplot as plt
def variance_fixed_NA():
    fs=open("./plot_data/fluctuation_of_entropy_Samplingtime100000RandomNA4.txt","r")
    data=fs.read().split()
    x=data[0::4]
    y=data[2::4]
    x=[float(i) for i in x]
    y=[float(i) for i in y]
    xmin=min(x)
    xmax=max(x)
    x2=np.linspace(xmin,xmax,1000)
    z2=-2.033*x2+4.983
    plt.figure(figsize=(12,8),dpi=300)
    plt.plot(x2,z2)
    plt.scatter(x,y,s=15,c='black')
    plt.xlim([3.5,5.5])
    plt.ylim([-6.5,-2])
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.legend("Fitting Curve",fontsize=40)
    plt.xlabel(r"$\ln N",fontsize=40)
    plt.ylabel(r"$\ln \sigma_{S_{A}}^{2}$",fontsize=40)
    plt.savefig("../script/figs/variance_fixed_NA4Period1.png")

variance_fixed_NA()