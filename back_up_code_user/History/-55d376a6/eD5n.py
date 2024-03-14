import matplotlib.pyplot as plt
import numpy as np

fig,ax=plt.subplots()
t=np.linspace(0,np.pi,100)
x=np.cos(t)
y=np.sin(t)*2/np.pi
ax.plot(x,y)
plt.xlabel(r'$w-v$')
plt.ylabel(r'$k_BT$')
plt.ylim(0,1)
ax.set_xticks([-1,0,1])
ax.set_xticklabels([r'$-u$','PT-broken',r'$u$'])
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')
ax.set_yticks([])
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
plt.show()