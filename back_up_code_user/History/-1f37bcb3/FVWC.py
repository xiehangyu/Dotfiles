import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
    N=int(sys.argv[1])
    temparray=np.zeros([N,N])
    temparray[0,0]=1
    for i in range(1,N):
        for j in range(N):
            startindex=min(j-1,0)
            temparray[i,j]=np.sum(temparray[i-1,startindex:])
    indexs=[]
    sums=[]
    logsums=[]
    for i in range(N):
        indexs.append(i+1)
        sums.append(np.sum(temparray[i,:]))
        logsums.append(np.log2(np.sum(temparray[i,:])))
    plt.plot(indexs,sums)
    plt.show()
    plt.plot(indexs,logsums)
    plt.show()