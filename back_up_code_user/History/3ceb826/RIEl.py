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
x1=[]
y1=[]
for i in range(len(xs)-1):
    x=np.linspace(xs[i],xs[i+1],100)
    y=a3[i]*x**3+a2[i]*x**2+a1[i]*x+a0[i]
    x1=x1+list(x)
    y1=y1+list(y)
    plt.plot(x,y,linewidth=3)
plt.plot(xs,ys,'o')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
ys[9]=10
gp2=open("./data/interpolationexpression_change.txt",'r')
gp2=gp2.read().split()
a3=[float(i) for i in gp2[4::4]]
a2=[float(i) for i in gp2[5::4]]
a1=[float(i) for i in gp2[6::4]]
a0=[float(i) for i in gp2[7::4]]
x2=[]
y2=[]
for i in range(len(xs)-1):
    x=np.linspace(xs[i],xs[i+1],100)
    y=a3[i]*x**3+a2[i]*x**2+a1[i]*x+a0[i]
    x2=x2+list(x)
    y2=y2+list(y)
    plt.plot(x,y,linewidth=3)
plt.plot(xs,ys,'o')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
plt.plot(x1,y1,'r',linewidth=3)
plt.plot(x2,y2,'b',linewidth=3)
plt.legend(['origin','change'])
plt.xlabel('x')
plt.ylabel('y')
plt.show()
