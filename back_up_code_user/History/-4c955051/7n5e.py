import numpy as np
import matplotlib.pyplot as plt
import scipy


def read_overlap_data(file_string):
    with open("./data/"+file_string+".csv",'r') as file:
        lines=file.readlines()
    n = np.min([len(lines),100])
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
    cutoff=1e-4
    below_curoff=np.sum(eigenvalues<cutoff)/len(eigenvalues)
    plt.hist(eigenvalues,bins=20,edgecolor='black')
    plt.xlabel("EigenValue")
    plt.ylabel("Frequency")
    plt.title("Rate={},cutoff={}".format(below_curoff,cutoff))
    plt.savefig("./figs/"+file_string+".png")
    plt.clf()
    return below_curoff


def complete_read_overlap_function(XX,N,alpha,t,tend,GSE,Korder,randomin,J):
    if XX==True:
        file_string="XXOBCOverlap"
    else:
        file_string="OBCOverlap"
    file_string=file_string+"N{}alpha{:.6f}t{:0.6f}tend{:.6f}GSE{}Korder{}randomin{:.6f}J{:.6f}".format(N,alpha,t,tend,GSE,Korder,randomin,J)
    return read_overlap_data(file_string)

if __name__=='__main__':
    t=0.03
    tend=9
    GSE=35
    Korder=5
    randomin=0.0
    for N in [40,60,80]:
        fp=open("./figs/XXOBCbelowcutoff{}N{}t{}GSE{}Korder{}randomin{}.txt".format(1e-4,N,t,GSE,Korder,randomin),'w')
        fp.write("Column alpha, row J\n")
        alphas=[0,0.1,0.2,0.4]
        fp.write("\t")
        for alpha in alphas:
            fp.write("{}\t".format(alpha))
        fp.write("\n")
        for alpha in alphas:
            for J in [0,0.2,0.5,0.8]:
                fp.write("{}\t".format(J))
                belowrate=complete_read_overlap_function(True,N,alpha,t,tend,GSE,Korder,randomin,J)
                fp.write("{}\t".format(belowrate))
            fp.write("\n")
    for N in [40,60,80]:
        fp=open("./figs/OBCbelowcutoff{}N{}t{}GSE{}Korder{}randomin{}.txt".format(1e-4,N,t,GSE,Korder,randomin),'w')
        fp.write("Column alpha, row J\n")
        alphas=[0,0.1,0.2,0.4]
        fp.write("\t")
        for alpha in alphas:
            fp.write("{}\t".format(alpha))
        fp.write("\n")
        for alpha in alphas:
            for J in [0,0.2,0.5,0.8]:
                fp.write("{}\t".format(J))
                belowrate=complete_read_overlap_function(True,N,alpha,t,tend,GSE,Korder,randomin,J)
                fp.write(belowrate+"\t")
            fp.write("\n")
