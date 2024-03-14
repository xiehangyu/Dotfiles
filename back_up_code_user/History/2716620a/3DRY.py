import numpy as np
import matplotlib.pyplot as plt

def draw_distance(N,alpha):
    fp=open("distancelongrange_N_"+str(N)+"_alpha_"+str(alpha)+".txt","r")
    data=fp.read().split()
    times=[float(i) for i in data[0::4]]
    


if __name__=='__main__':
    decision=input("1 for distance and drawing\n 2 for ratio drawing\n")
    N=input("input N\n")
    alpha=input("input alpha\n")
    if decision==1:
        draw_distance(N,alpha)