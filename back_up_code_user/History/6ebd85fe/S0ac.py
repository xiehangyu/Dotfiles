import numpy as np
import matplotlib.pyplot as plt
x=[32,16,8]
y1=[14235069,19185059,29589045]
y2=[15104757,15471124,22050610]
y1=np.array(y1)
y2=np.array(y2)
y1=np.log2(y1)
y2=np.log2(y2)
plt.rcParams.update({'font.size': 22})
plt.figure(figsize=(12,8))
plt.plot(x,y1,'r')
plt.plot(x,y2,'b')
plt.legend(['Basic','Shared Memory','',''])
plt.scatter(x,y1,c='r')
plt.scatter(x,y2,c='b')
plt.xlabel(r"Block Size")
plt.ylabel(r'$\log_2$(Time)')
plt.title(r"N=2048")
plt.tight_layout()
plt.savefig('./Figs/GPU.png',dpi=300)
plt.show()
