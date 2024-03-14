import numpy as np
import scipy
from floquet_hamiltonian import foureir_transform_from_kspace_to_real_space_period1

class dynamical_monte:
    def __init__(self,N,norder,SamplingSize):
        self.Nsize=N
        self.norder=norder
        self.Samplingsize=SamplingSize
        self.fourier_transformation_from_kspace_to_real_space=foureir_transform_from_kspace_to_real_space_period1(N)
        indexs=[i for i in range(N)]
        indexs=np.array(indexs)
        momentums=2*np.pi/N*indexs
        self.energies=2*np.cos(momentums)
        self.init_covariance_matrix_position=np.zeros([N,N])
        for i in range(0,N,2):
            self.init_covariance_matrix_position[i,i]=1
        self.init_covariance_matrix_momentum=self.fourier_transformation_from_kspace_to_real_space.T.conj().dot(self.init_covariance_matrix_position.dot(self.fourier_transformation_from_kspace_to_real_space))       