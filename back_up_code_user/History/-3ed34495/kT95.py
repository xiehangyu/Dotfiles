import numpy as np
import matplotlib.pyplot as plt

fp=open("./data/iris_PCA.txt","r")
data=fp.read().split()
x=[float(i) for i in data[0::3]]
y=[float(i) for i in data[1::3]]
label=[int(i) for i in data[2::3]]
colors = np.array(['r', 'g', 'b'])
label_colors = colors[label]
fig,ax=plt.subplots()
scatter=ax.scatter(x,y,c=label_colors)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='0 - Red',
                          markerfacecolor='r', markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='1 - Green',
                          markerfacecolor='g', markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='2 - Blue',
                          markerfacecolor='b', markersize=10)]
ax.legend(handles=legend_elements, loc='lower left')

plt.show()
plt.show()