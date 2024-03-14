import numpy as np
import scipy
from floquet_hamiltonian import foureir_transform_from_kspace_to_real_space_period1

class dynamical_monte:
    def __init__(self,N,norder,SamplingSize):
        self.Nsize=N
        self.norder=norder
        self.Samplingsize=SamplingSize
        self.fourier_transformation_from_kspace_to_real_space
        