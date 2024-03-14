from floquet_hamiltonian import foureir_transform_from_kspace_to_real_space as fourier_transform_from_kspace_to_real_space_period2
from scipy.linalg import block_diag
from floquet_hamiltonian import generate_covariance_matrix_density_wave
from floquet_hamiltonian import generate_covariance_matrix_half
import numpy as np

def generate_Hamiltonian_in_AB_basis(N_site):
    Hamiltonian_list=[]
    J1=0.3
    J2=-0.3
    for i in range(0,N_site//2):
        k=2*np.pi*i/N_site
        temp_matrix=np.array([[0.5,2*np.cos(k)],[2*np.cos(k),0]])
        #temp_matrix=np.array([[0,J1*np.exp(5j*k)+J2*np.exp(-5j*k)+2*np.cos(k)],[J1*np.exp(-5j*k)+J2*np.exp(5j*k)+2*np.cos(k),0]])
        Hamiltonian_list.append(temp_matrix)
    return block_diag(*Hamiltonian_list)


def calculate_covariance_matrix_in_conserved_basis(N_site):
    covariance_matrix=generate_covariance_matrix_half(Nsite=N_site)
    print(np.diag(covariance_matrix))
    UF=fourier_transform_from_kspace_to_real_space_period2(N_site)
    covariance_matrix_in_AB_basis=UF.T.conj().dot(covariance_matrix.dot(UF))
    Hamiltonian_in_AB_basis=generate_Hamiltonian_in_AB_basis(N_site) 
    v,u=np.linalg.eig(Hamiltonian_in_AB_basis)
    return u.T.dot(covariance_matrix_in_AB_basis.dot(u.conj()))
    
def transform_matrix_from_real_basis_to_conserved_basis(N_site):
    UF=fourier_transform_from_kspace_to_real_space_period2(N_site)
    Hamiltonian_in_AB_basis=generate_Hamiltonian_in_AB_basis(N_site)
    v,u=np.linalg.eig(Hamiltonian_in_AB_basis)
    return u.T.dot(UF.T.conj())


def Q_conserve(i,N_site):
    temp_matrix=transform_matrix_from_real_basis_to_conserved_basis(N_site)
    temp_array=temp_matrix[2*i,:]
    temp_array=np.array([list(temp_array)])
    return temp_array.T.conj().dot(temp_array)

def P_conserve(i,N_site):
    temp_matrix=transform_matrix_from_real_basis_to_conserved_basis(N_site)
    temp_array=temp_matrix[2*i+1,:]
    temp_array=np.array([list(temp_array)])
    return temp_array.T.conj().dot(temp_array)