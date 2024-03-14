import matplotlib.pyplot as plt
import numpy as np

t=np.linspace(0,np.pi,100)
x=np.sin(t)
y=np.cos(t)
plt.plot(x,y)
plt.show()