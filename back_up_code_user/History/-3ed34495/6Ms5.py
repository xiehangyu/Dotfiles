import numpy as np
import matplotlib.pyplot as plt

fp=open("./data/iris_PCA.txt","r")
data=fp.read().split()
x=[float(i) for i in data[0::3]]
y=[float(i) for i in data[1::3]]
label=[int(i) for i in data[2::3]]
colors = np.array(['r', 'g', 'b'])
label_colors = colors[label]
scatter=plt.scatter(x,y,c=label_colors)
plt.legend(*scatter.legend_elements(),title="Label")
plt.show()