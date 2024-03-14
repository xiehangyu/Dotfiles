import matplotlib.pyplot as plt
import numpy as np
import matplotlib.axes

fig,ax=plt.subplots()
t=np.linspace(0,np.pi,100)
x=np.cos(t)
y=np.sin(t)*2/np.pi
ax.plot(x,y)
plt.xlabel(r'$w-v$',fontsize=20)
plt.ylabel(r'$k_BT$',fontsize=20)
plt.ylim(0,1)
ax.set_xticks([-1,0,1])
ax.set_xticklabels([r'$-u$','PT-broken',r'$u$'])
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')
ax.set_yticks([])
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(1.05, 0.1)
ax.yaxis.set_label_coords(0.60, 0.975)
xmin, xmax = ax.get_xlim()
ax.arrow(xmax, 0, 0.2, 0.0, head_width = 0.05, 
         head_length = 0.1)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()