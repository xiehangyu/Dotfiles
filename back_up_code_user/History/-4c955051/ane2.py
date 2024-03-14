import numpy as np
import matplotlib.pyplot as plt


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