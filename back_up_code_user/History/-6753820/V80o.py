import numpy as np
from scipy.optimize import curve_fit

def my_func1(x,a1,a2,a3,a4):
    return a1*(1-np.exp(-x^a2/a3))+a4

def my_func2(x,a1,a2,a3,a4,a5):
    return a1*np.exp(-x/a2)+a3*np.exp(-x/a4)+a5

