import numpy as np
import matplotlib.pyplot as plt


def read_overlap_data(file_string):
    with open("./data/"+file_string+".csv",'r') as file:
        lines=file.readlines()