import numpy as np 
import matplotlib.pyplot as plt

def minimum_distance(i,j,N):
    return min(abs(i-j),N-abs(i-j))
def construct_matrix(N,lamb):
    matrix=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if i==j:
                matrix[i,j]=-1
            else:
                matrix[i,j]=-lamb**minimum_distance(i,j,N)
    eigenvalues=np.sort(np.linalg.eigvals(matrix))
    eigenvalues=eigenvalues-eigenvalues[0]
    number=np.min([N,20])
    return eigenvalues[:number]

def draw_plot(N):
    eigs=[]
    lambs=np.linspace(0,1,40)
    for lamb in lambs:
        eigs.append(construct_matrix(N,lamb))
    eigs=np.array(eigs)
    for i in range(eigs.shape[1]):
        plt.plot(lambs,eigs[:,i])
        plt.scatter(lambs,eigs[:,i])
    plt.show()

if __name__=="__main__":
    draw_plot(50)

