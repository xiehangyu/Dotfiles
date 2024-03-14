import numpy as np
import matplotlib.pyplot as plt
import scipy


def read_overlap_data(file_string):
    with open("./data/"+file_string+".csv",'r') as file:
        lines=file.readlines()
    n = len(lines)
    matrix=np.zeros((n,n),dtype=complex)
    row=0
    for line in lines:
        entries=line.split(',')
        complex_entries=[complex(float(entries[i]),float(entries[i+1])) for i in range(0,len(entries),2)]
        matrix[row,:len(complex_entries)]=complex_entries
        row+=1
    for i in range(n):
        for j in range(i+1,n):
            matrix[i,j]=np.conj(matrix[j,i])
    eigenvalues=scipy.linalg.eigvalsh(matrix)
    bins=np.arange(0,1.1,0.1)
    plt.hist(eigenvalues,bins=bins,edgecolor='black')
    plt.xlabel("EigenValue")
    plt.ylabel("Frequency")
    plt.savefig("./figs/"+file_string+".png")


def complete_read_overlap_function(XX,N,alpha,t,tend,GSE,Korder,randomin,J):
    if XX==True:
        file_string="XXOBCOverlap"
    else:
        file_string="OBCOverlap"
    file_string=file_string+"N{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}.csv".format(N,alpha,t,tend,GSE,Korder,randomin,J)
    read_overlap_data(file_string)

if __name__=='__main__':
    t=0.1
    tend=30
    GSE=35
    Korder=5
    randomin=0.0
    N=80
    for alpha in [0,0.1]:
        for J in [0,0.2,0.5,0.8]:
            complete_read_overlap_function(True,N,alpha,t,tend,GSE,Korder,randomin,J)
    
