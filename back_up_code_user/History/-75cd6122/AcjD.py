import numpy as np
import matplotlib.pyplot as plt
def variance_fixed_NA():
    fs=open("./plot_data/fluctuation_of_entropy_Samplingtime100000RandomNA4.txt","r")
    data=fs.read().split()
    x=data[0::4]
    y=data[2::4]
    x=[float(i) for i in x]
    y=[float(i) for i in y]

