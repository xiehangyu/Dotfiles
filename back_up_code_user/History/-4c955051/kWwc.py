import numpy as np
import matplotlib.pyplot as plt


def read_overlap_data(file_string):
    with open("./data/"+file_string+".csv",'r') as file:
        lines=file.readlines()
    n = len(lines)
    matrix=np.zeros((n,n),dtype=complex)
    for line in lines:
        entries=line.split(',')
        complex_entries=[complex(float(entries[i],float(entries[i+1]))) for i in range(0,len(entries),2)]