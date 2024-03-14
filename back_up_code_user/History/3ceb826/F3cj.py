import numpy as np
import matplotlib.pyplot as plt
fp1=open("./data/point.txt",'r')
gp1=open("./data/interpolationexpression_origin.txt",'r')
fp1=fp1.read().split()
xs=[float(i) for i in fp1[0::2]]
ys=[float(i) for i in fp1[1::2]]
gp1=gp1.read().split()
a3=[float(i) for i in gp1[4::4]]
a2=[float(i) for i in gp1[5::4]]
a1=[float(i) for i in gp1[6::4]]
a0=[float(i) for i in gp1[7::4]]
for i in range(len(xs)-1):
    x=np.linspace(xs[i],xs[i+1],100)
    y=a3[i]*x**3+a2[i]*x**2+a1[i]*x+a0[i]
    plt.plot(x,y,linewidth=3)
plt.plot(xs,ys,'o')
plt.show()