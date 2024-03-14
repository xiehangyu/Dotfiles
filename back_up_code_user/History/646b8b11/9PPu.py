import numpy as np
import matplotlib.pyplot as plt
fp=open("./data/off_diagonal_square1.txt","r")
data=fp.read().split()
x=[float(i) for i in data[0::2]]
y=[float(i) for i in data[1::2]]
plt.plot(x,y)
plt.scatter(x,y)
plt.xlabel('step',fontsize=30)
plt.ylabel('off-diagonal-square',fontsize=30)
plt.tick_params(axis='x', labelsize=30)
plt.tick_params(axis='y', labelsize=30)
plt.show()