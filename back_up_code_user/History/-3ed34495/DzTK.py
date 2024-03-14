fp=open("./data/iris_PCA.txt","r")
data=fp.read().split()
x=[float(i) for i in data[0::3]]
y=[float(i) for i in data[1::3]]
label=[int(i) for i in data[2::3]]
colors=["red","green","blue"]
