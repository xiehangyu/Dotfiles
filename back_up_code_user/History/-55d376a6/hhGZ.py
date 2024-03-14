import matplotlib.pyplot as plt
import numpy as np

fig,ax=plt.subplots()
t=np.linspace(0,np.pi,100)
x=np.cos(t)
y=np.sin(t)*2/np.pi
ax.plot(x,y)
plt.xlabel(r'$w-v$')
plt.ylabel(r'$k_BT$')
ax.set_xticks([-1,0,1])
ax.set_xticklabels([r'$-u$','PT-broken',r'$u$'])
plt.show()