from first_trial import *
from matplotlib import cm
from scipy.linalg import expm
from random_hamiltonian import *
from multiprocessing import Process, Array
from measure_concentration import hilbert_schmidt_norm
import sys
import os
import warnings
import time
import scipy.stats
from scipy.linalg import block_diag
from random import randint
import random
import datetime

warnings.filterwarnings('error')
gauss_sigma=1/3
N_sites=500
T_period=0.0004
T_totaltime=100000000
Sampling_time=250000
GLOBAL_METHOD=1
sigma_x=np.array([[0,1],[1,0]])
sigma_y=np.array([[0,-1j],[1j,0]])
sigma_z=np.array([[1,0],[0,-1]])
decay_length=25

def Hamiltonian1():
    if GLOBAL_METHOD==2:
        return Hamiltonian1_disuniform()
    H1=np.zeros((N_sites,N_sites))
    for i in range(0,N_sites,2):
        H1[i,(i+1)%N_sites]=1
        H1[(i+1)%N_sites,i]=1
    return H1

def Hamiltonian1_disuniform():
    H1=np.zeros((N_sites,N_sites))
    for i in range(0,N_sites,2):
        temp=np.random.normal(1,gauss_sigma)
        H1[i,(i+1)%N_sites]=temp
        H1[(i+1)%N_sites,i]=temp
    return H1


def Hamiltonian2():
    if GLOBAL_METHOD==2:
        return Hamiltonian2_disuniform()
    H2=np.zeros((N_sites,N_sites))
    for i in range(1,N_sites,2):
        H2[i,(i+1)%N_sites]=1
        H2[(i+1)%N_sites,i]=1
    return H2

def Hamiltonian3():
    H3=np.zeros((N_sites,N_sites),dtype=complex)
    for i in range(N_sites):
        phase=np.random.uniform(0,2*np.pi)
        H3[i,(i+1)%N_sites]=1
        H3[(i+1)%N_sites,i]=1
        if i%2==0:
            H3[i,(i+2)%N_sites]+=0.0
            H3[(i+2)%N_sites,i]+=0.0
            #H3[i,(i+3)%N_sites]=0.0
            #H3[(i+3)%N_sites,i]=0.0
        if i%2==1:
            H3[i,(i+2)%N_sites]+=+0.0
            H3[(i+2)%N_sites,i]+=+0.0
            #H3[i,(i+3)%N_sites]=-0.0
            #H3[(i+3)%N_sites,i]=-0.0
#    H3=np.load("./plot_data/M_matrix_Tperiod{}Nsites{}.npy".format(T_period,N_sites))
    return H3


def Hamiltonian_decay():
    H4=np.zeros((N_sites,N_sites),dtype=complex)
    for i in range(N_sites):
        for j in range(i+1,N_sites):
            distance=min(abs(i-j),N_sites-abs(i-j))
            value=np.exp(-(distance-1)/decay_length)
            phase=np.random.uniform(0,2*np.pi)
            H4[i,j]=value*np.exp(phase*1j)
            H4[j,i]=value*np.exp(-phase*1j)
    return H4

def Hamiltonian2_disuniform():
    H2=np.zeros((N_sites,N_sites))
    for i in range(1,N_sites,2):
        temp=np.random.normal(1,gauss_sigma)
        H2[i,(i+1)%N_sites]=temp
        H2[(i+1)%N_sites,i]=temp
    return H2

def evolution(original_matrix,evolution_matrix):
    after_matrix=expm(1j*T_period*evolution_matrix).dot(original_matrix)
    return after_matrix


def entropy_sampling(covariance_matrix, starting_point, sampling_sites):
    calculating_matrix=covariance_matrix[starting_point:starting_point+sampling_sites,starting_point:starting_point+sampling_sites]
    return calculate_entropy_from_correlation_matrix(calculating_matrix)


def page_curve_over_time_average(calculate_list):
    covariance_matrix=np.zeros((N_sites,N_sites))
    H1=Hamiltonian1()
    H2=Hamiltonian2()
    H3=Hamiltonian3()
    H4=Hamiltonian_decay()
    for i in range(0,N_sites,2):
        covariance_matrix[i,i]=1
    A_matrix=np.eye(N_sites)
    for i in range(T_totaltime):
        if i%100==0:
            print(i)
       # A_matrix=evolution(A_matrix,H1)
       # A_matrix=evolution(A_matrix,H2)
        A_matrix=evolution(A_matrix, H3)
    entropy_array=np.zeros(len(calculate_list))
    for j in range(Sampling_time):
        if j%100==0:
            print(j)
        later_time_covariance_matrix=A_matrix.T.dot(covariance_matrix.dot(A_matrix.conj()))
        for k in range(len(calculate_list)):
            entropy_array[k]=entropy_array[k]+entropy_sampling(later_time_covariance_matrix,N_sites//2-calculate_list[k]//2,calculate_list[k])
        #A_matrix=evolution(A_matrix,H1)
        #A_matrix=evolution(A_matrix,H2)
        A_matrix=evolution(A_matrix,H3)
    return entropy_array/Sampling_time

def page_curve_comparison():
    calculate_list=[i for i in range(1,N_sites,2)]
    entropy_array=page_curve_over_time_average(calculate_list)
    entropy_array2=[]
    fs=open("./plot_data/Nsite{}_Tperiod{}_Totaltime{}_Method{}Samplingtime{}.txt".format(N_sites,T_period,T_totaltime,GLOBAL_METHOD,Sampling_time),'w')
    for j in range(len(calculate_list)):
        fs.write("{}\t{}\n".format(calculate_list[j],entropy_array[j]))
        instance1=trun_k_probability(N_sites,N_sites//2 , calculate_list[j])
        entropy_array2.append(instance1.direct_entropy_without_sum())
    fs.close()
    plt.plot(calculate_list,entropy_array,calculate_list,entropy_array2)
    plt.legend(["Local Floquet Hamiltonian","Random State"])
    #plt.show()

def periodic_of_field_string():
    global T_period
    calculate_list=[5]
    T_periods=[]
    entropys=[]
    for i in list(np.linspace(np.pi/2,np.pi,20)):
        T_periods.append(i)
        T_period=i
        print(T_period)
        entropys.append(page_curve_over_time_average(calculate_list)[0])
    fs=open("./plot_data/periodic_of_field_strength_method{}.txt".format(GLOBAL_METHOD),"w")
    for j in range(len(T_periods)):
        fs.write("{}\t{}\n".format(T_periods[j],entropys[j]))
    fs.close()
    plt.plot(T_periods,entropys)
    plt.show()

def fluctuation_of_entropy(initial_time):
    global N_sites
    N_input=sys.argv[1]
    N_input=int(N_input)
    N_sites=N_input
    H3=Hamiltonian3()
    H1=Hamiltonian1()
    H2=Hamiltonian2()
    covariance_matrix=generate_covariance_matrix_density_wave(Nsite=N_input)
    print(np.diag(covariance_matrix))
    A_matrix=np.eye(N_sites)
    UF=foureir_transform_from_kspace_to_real_space(Nsite=N_input)
    #A_matrix=evolution_matrix_in_real_space(initial_time, UF).T
    A_matrix=expm(1j*initial_time*H3).dot(A_matrix)
    entropy_square_sum=0
    entropy_sum=0
    min_value=1000
    max_value=0
    for i in range(Sampling_time):
        calculate_matrix=A_matrix.T.dot(covariance_matrix.dot(A_matrix.conj()))
        A_matrix=evolution(A_matrix,H3)
    #    A_matrix=evolution(A_matrix,H2)
        entropy=entropy_sampling(calculate_matrix,N_sites//2-2,4)
        entropy_sum+=entropy
        entropy_square_sum+=entropy*entropy
        if(max_value<entropy):
            max_value=entropy
        if(min_value>entropy):
            min_value=entropy
    entropy_sum=entropy_sum/Sampling_time
    entropy_square_sum=entropy_square_sum/Sampling_time
    entropy_square_sum=entropy_square_sum-entropy_sum**2
    fs=open("./plot_data/fluctuation_of_entropy_Samplingtime{}period2Range3NA4.txt".format(Sampling_time),"a+")
    fs.write("{}\t{}\t{}\t{}\n".format(N_sites,entropy_sum,entropy_square_sum,max_value-min_value))
    fs.close()


def read_equilibrium(initial_time_step=0,mask_matrix=None):
    covariance_matrix=np.zeros((N_sites,N_sites))
    H1=Hamiltonian1()
    H2=Hamiltonian2()
    H3=Hamiltonian3()
    H4=Hamiltonian_decay()
    covariance_matrix=generate_covariance_matrix_density_wave()
    print(np.diag(covariance_matrix))
    #covariance_matrix=np.eye(N_sites)/2
    if initial_time_step==0:
        A_matrix=np.eye(N_sites)
    else:
        A_matrix=np.eye(N_sites)
        A_matrix=expm(1j*initial_time_step*H3).dot(A_matrix)
        UF=foureir_transform_from_kspace_to_real_space()
        #A_matrix=evolution_matrix_in_real_space(initial_time_step, UF).T
    times=[]
    entropys=[]
    i=0
    ###Define the reference line
    references=[]
    ###Define the reference line
    while True:
        i=i+1
        calculate_matrix=A_matrix.T.dot(covariance_matrix.dot(A_matrix.conj()))
    #    A_matrix=evolution(A_matrix,H1)
    #    A_matrix=evolution(A_matrix,H2)
        A_matrix=evolution(A_matrix,H3)
        if type(mask_matrix)!=type(None):
            temp_matrix=mask_matrix*calculate_matrix
            entropy=temp_matrix.sum()
            entropy=abs(entropy)
        else:
            entropy=entropy_sampling(calculate_matrix,N_sites//2-30,60)
        times.append(i*T_period)
        entropys.append(entropy)
        ###Draw the reference line
        sum_entropy_reference=0
        for iter_index in range(N_sites):
            velocity=2*np.sin(2*np.pi*iter_index/N_sites)
            velocity=np.abs(velocity)
            sum_entropy_reference += 2*i*T_period*velocity/N_sites
        references.append(sum_entropy_reference)
        ###Draw the reference line
        if i%100==0:
            plt.cla()
            plt.plot(times,entropys)
        ###Draw the reference line
            plt.plot(times,references)
            plt.legend(["Entropy","Reference"])
        ###Draw the reference line
            plt.pause(0.5)
            

def generate_floquet_E_n_in_kspace(k):
    if k==np.pi/2:
        return 0,0,0,0
    E=np.arccos(np.cos(T_period)**2-np.sin(T_period)**2*np.cos(2*k))/(2*T_period)
    try:
        nx=2*np.sin(T_period)*np.cos(T_period)*np.cos(k)/np.sin(2*T_period*E)
        ny=0
        nz=-np.sin(T_period)**2*np.sin(2*k)/np.sin(2*T_period*E)
    except Warning:
        print("{}\t{}\t{}\n".format(k,E,T_period))
        nx=0
        ny=0
        nz=0
    return E,nx,ny,nz


def generate_A_matrix_in_kspace(k,n_step):
    E,nx,ny,nz=generate_floquet_E_n_in_kspace(k)
    A_matrix=np.cos(2*n_step*T_period*E)*np.eye(2)+1j*np.sin(2*n_step*T_period*E)*(nx*sigma_x+ny*sigma_y+nz*sigma_z)
    return A_matrix


def generate_M_matrix_in_kspace(k):
    E,nx,ny,nz=generate_floquet_E_n_in_kspace(k)
    M_matrix=E*(nx*sigma_x+ny*sigma_y+nz*sigma_z)
    return M_matrix


def generate_full_A_matrix_in_kspace(n_step):
    block_elements=[]
    for i in range(N_sites//2):
        block_elements.append(generate_A_matrix_in_kspace(2*np.pi*i/N_sites,n_step))
    return block_diag(*block_elements)

def generate_full_M_matrix_in_kspace():
    block_elements=[]
    for i in range(N_sites//2):
        block_elements.append(generate_M_matrix_in_kspace(2*np.pi*i/N_sites))
    return block_diag(*block_elements)



def foureir_transform_from_kspace_to_real_space(Nsite=N_sites):
    UF=np.zeros((Nsite,Nsite),dtype=complex)
    for i in range(Nsite):
        k=2*np.pi/Nsite*(i//2)
        if i%2==0:
            for l in range(0,Nsite,2):
                UF[i,l]=np.exp(-1j*k*(l+1))
        if i%2==1:
            for l in range(1,Nsite,2):
                UF[i,l]=np.exp(-1j*k*(l+1))
    return np.sqrt(2/Nsite)*UF.T.conj()


def foureir_transform_from_kspace_to_real_space_period1(Nsites=N_sites):
    UF=np.zeros((Nsites,Nsites),dtype=complex)
    for i in range(Nsites):
        k=2*np.pi*i/Nsites
        for l in range(Nsites):
            UF[i,l]=np.exp(-1j*k*(l+1))
    return np.sqrt(1/Nsites)*UF.T.conj()

def evolution_matrix_in_real_space(n_step,UF):
    A_mat=generate_full_A_matrix_in_kspace(n_step)
    return UF.dot(A_mat.T.dot(UF.T.conj()))


def M_matrix_in_real_space(UF):
    M_mat=generate_full_M_matrix_in_kspace()
    return UF.conj().dot(M_mat.dot(UF.T))

def locality_of_evolution_matrix(step):
    UF=foureir_transform_from_kspace_to_real_space()
    Z=evolution_matrix_in_real_space(step, UF)
    Z=np.abs(Z)
    x=[i for i in range(N_sites)]
    y=[i for i in range(N_sites)]
    X,Y=np.meshgrid(x,y)
    plt.contourf(X,Y,Z,cmap='hot_r')
    plt.colorbar()
    plt.show()

def locality_of_Hamiltonian_decay():
    Z=Hamiltonian_decay()
    Z=np.abs(Z)
    x=[i for i in range(N_sites)]
    y=[i for i in range(N_sites)]
    X,Y=np.meshgrid(x,y)
    plt.contourf(X,Y,Z,cmap='hot_r')
    plt.colorbar()
    plt.show()

def locality_of_M_matrix():
    UF=foureir_transform_from_kspace_to_real_space()
    Z=M_matrix_in_real_space(UF)
    np.save("./plot_data/M_matrix_Tperiod{}Nsites{}".format(T_period,N_sites),Z,allow_pickle=False)
    Z=(np.real(Z))
    x=[i for i in range(N_sites)]
    y=[i for i in range(N_sites)]
    X,Y=np.meshgrid(x,y)
    plt.contourf(X,Y,Z,cmap='hot_r')
    plt.colorbar()
    plt.show()

def generate_covariance_matrix_random(space=2,Nsite=N_sites):
    random.seed(time.time())
    covariance_matrix=np.zeros([Nsite,Nsite])
    i=0
    while i<Nsite//space:
        i=i+1
        temp_index=randint(0,Nsite-1)
        if(covariance_matrix[temp_index,temp_index]==1):
            i=i-1
            continue
        else:
            covariance_matrix[temp_index,temp_index]=1
    #print(np.diag(covariance_matrix))
    return covariance_matrix

def generate_covariance_matrix_density_wave(space=2,Nsite=N_sites):
    covariance_matrix=np.zeros([Nsite,Nsite])
    for i in range(0,Nsite,space):
        covariance_matrix[i,i]=1
    #print(np.diag(covariance_matrix))
    return covariance_matrix 


def generate_covariance_matrix_half(space=2,Nsite=N_sites):
    covariance_matrix=np.zeros([Nsite,Nsite])
    for i in range(0,Nsite//space):
        covariance_matrix[i,i]=1
    #print(np.diag(covariance_matrix))
    return covariance_matrix

def page_curve_over_time_average2(calculate_list):
    UF=foureir_transform_from_kspace_to_real_space()
    covariance_matrix=np.zeros([N_sites,N_sites])
    covariance_matrix=generate_covariance_matrix_density_wave()
    print(np.diag(covariance_matrix))
    #i=0
    #while i<N_sites//2:
    #    i=i+1
    #    temp_index=randint(0,N_sites-1)
    #    if(covariance_matrix[temp_index,temp_index]==1):
    #        i=i-1
    #        continue
    #    else:
    #        covariance_matrix[temp_index,temp_index]=1
    entropy_array=np.zeros(len(calculate_list))
    time_clock=0
    for step in range(T_totaltime,T_totaltime+Sampling_time):
        time_clock+=1
        if (time_clock%100==0):
            print(time_clock)
        evolution_matrix=evolution_matrix_in_real_space(step,UF)
        after_matrix=evolution_matrix.dot(covariance_matrix.dot(evolution_matrix.T.conj()))
        for k in range(len(calculate_list)):
            entropy_array[k]=entropy_array[k]+entropy_sampling(after_matrix,N_sites//2-calculate_list[k]//2, calculate_list[k])    
    return entropy_array/Sampling_time
    

def page_curve_comparison2():
    calculate_list=[i for i in range(1,N_sites,2)]
    entropy_array=page_curve_over_time_average2(calculate_list)
    entropy_array2=[]
    fs=open("./plot_data/Nsite{}_Tperiod{}_Totaltime{}_Method{}MF.txt".format(N_sites,T_period,T_totaltime,GLOBAL_METHOD),'w')
    for j in range(len(calculate_list)):
        fs.write("{}\t{}\n".format(calculate_list[j],entropy_array[j]))
        instance1=trun_k_probability(N_sites,N_sites//2 , calculate_list[j])
        entropy_array2.append(instance1.direct_entropy_without_sum())
    fs.close()
    plt.plot(calculate_list,entropy_array,calculate_list,entropy_array2)
    plt.legend(["Local Floquet Hamiltonian","Random State"])
    #plt.show()
    
def measure_different_T_period_locality(c1):
    global T_period
    xs=[i for i in range(N_sites)]
    yss=[]
    UF=foureir_transform_from_kspace_to_real_space()
    for T_period in c1:
        Z=M_matrix_in_real_space(UF)
        Z=np.abs(Z)
        yss.append(Z[N_sites//2])
    for ys in yss:
        plt.plot(xs,ys)
    plt.legend([str(char) for char in c1])
    plt.show()

def hilbert_schmidt_norm_with_T_period():
    global T_period
    UF=foureir_transform_from_kspace_to_real_space()
    xs=[]
    ys=[]
    for T_period in np.linspace(0,np.pi,100):
        if T_period ==0:
            continue
        Z=M_matrix_in_real_space(UF)
        xs.append(T_period)
        ys.append(hilbert_schmidt_norm(Z))
    fs=open('./plot_data/norm_of_M_with_T_period_Nsites{}.txt'.format(N_sites),'w')
    for j in range(len(xs)):
        fs.write("{}\t{}\n".format(xs[j],ys[j]))
    plt.plot(xs,ys)
    plt.show()


def page_curve_over_time_average_hamiltonian_evolution(calculate_list,initial_time,spacing):
    covariance_matrix=np.zeros((N_sites,N_sites))
    H3=Hamiltonian3()
    covariance_matrix=generate_covariance_matrix_density_wave(spacing)
    print(np.diag(covariance_matrix))
    A_matrix=np.eye(N_sites)
    A_matrix=expm(1j*initial_time*H3).dot(A_matrix)
    entropy_array=np.zeros(len(calculate_list))
    for j in range(Sampling_time):
        if j%100==0:
            print(j)
        later_time_covariance_matrix=A_matrix.T.dot(covariance_matrix.dot(A_matrix.conj()))
        for k in range(len(calculate_list)):
            entropy_array[k]=entropy_array[k]+entropy_sampling(later_time_covariance_matrix,N_sites//2-calculate_list[k]//2,calculate_list[k])
        #A_matrix=evolution(A_matrix,H1)
        #A_matrix=evolution(A_matrix,H2)
        A_matrix=evolution(A_matrix,H3)
    return entropy_array/Sampling_time




def page_curve_comparison_hamiltonian_evolution(initial_time,spacing=2):
    calculate_list=[i for i in range(1,N_sites,2)]
    entropy_array=page_curve_over_time_average_hamiltonian_evolution(calculate_list,initial_time,spacing)
    entropy_array2=[]
    fs=open("./plot_data/Hamiltonian_evolution_Nsite{}_Tperiod{}_Initialtime{}_Samplingtime{}Spacing{}Period1densitywave.txt".format(N_sites,T_period,initial_time,Sampling_time,spacing),'w')
    for j in range(len(calculate_list)):
        fs.write("{}\t{}\n".format(calculate_list[j],entropy_array[j]))
        instance1=trun_k_probability(N_sites,N_sites//spacing , calculate_list[j])
        entropy_array2.append(instance1.direct_entropy_without_sum())
    fs.close()
    plt.plot(calculate_list,entropy_array,calculate_list,entropy_array2)
    plt.legend(["Local Floquet Hamiltonian","Random State"])
    plt.show()


def conserved_quantity_of_floquet_system():
    UF=foureir_transform_from_kspace_to_real_space()
    Z=M_matrix_in_real_space(UF)
    w,v=np.linalg.eig(Z)
    return v.T.conj()

def Random_Gaussian_state_variance():
    Nsite=sys.argv[1]
    Nsite=int(Nsite)
    covariance_matrix=generate_covariance_matrix_density_wave(space=2,Nsite=Nsite)
    entropy_sum=0
    entropy_square_sum=0
    min_value=1000
    max_value=0
    for i in range(Sampling_time):
        unitary=scipy.stats.unitary_group.rvs(Nsite)
        calculate_matrix=unitary.dot(covariance_matrix.dot(unitary.T.conj()))
        entropy=entropy_sampling(calculate_matrix, Nsite//2-2, 4)
        entropy_sum += entropy
        entropy_square_sum += entropy*entropy
        if entropy<min_value:
            min_value=entropy
        if entropy>max_value:
            max_value=entropy
    entropy_sum=entropy_sum/Sampling_time
    entropy_square_sum=entropy_square_sum/Sampling_time
    entropy_square_sum=entropy_square_sum-entropy_sum*entropy_sum
    fs=open("./plot_data/fluctuation_of_entropy_Samplingtime{}RandomNA4.txt".format(Sampling_time),"a+")
    fs.write("{}\t{}\t{}\t{}\n".format(Nsite,entropy_sum,entropy_square_sum,max_value-min_value))
    fs.close()

if __name__=='__main__':
    #page_curve_comparison2()
    mask_matrix=np.zeros([N_sites,N_sites])
    for i in range(N_sites):
        mask_matrix[i,(i+1)%N_sites]=1
        mask_matrix[(i+1)%N_sites,i]=1
    #read_equilibrium(0,mask_matrix=None)
    page_curve_comparison_hamiltonian_evolution(100000000,2)
    #page_curve_comparison2()
    #Random_Gaussian_state_variance()
    #fluctuation_of_entropy(100000000)
