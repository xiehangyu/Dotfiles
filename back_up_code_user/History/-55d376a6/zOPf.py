import matplotlib.pyplot as plt
import numpy as np

t=np.linspace(0,np.pi,100)
x=np.cos(t)
y=np.sin(t)*2/np.pi
plt.plot(x,y)
plt.show()
plt.xlabel(r'$w-v$')
plt.ylabel(r'$k_BT$')
plt.