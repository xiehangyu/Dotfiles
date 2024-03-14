from time import time
from first_trial import *
from local_conservable_analytical_calculation import transform_matrix_from_real_basis_to_conserved_basis
from floquet_hamiltonian import generate_covariance_matrix_density_wave
def read_page_curve_comparison():
    N_sites=input("N_sites:")
    N_sites=int(N_sites)
    Totaltime=input("Totaltime:")
    Samplingtime=input("Sampling Time:")
    T_period=input("T_period:")
    fs=open("./plot_data/Nsite{}_Tperiod{}Totaltime{}Samplingtime{}Create_byC.txt".format(N_sites,T_period,Totaltime,Samplingtime),"r")
    data=fs.read().split()
    fs.close()
    Calculate_list=data[0::2]
    value=data[1::2]
    Calculate_list=[int(i) for i in Calculate_list]
    value=[float(i) for i in value]
    entropy_array2=[]
    for j in Calculate_list:
        instance1=trun_k_probability(N_sites,N_sites//2,j)
        entropy_array2.append(instance1.direct_entropy_without_sum())
    plt.plot(Calculate_list,value,Calculate_list,entropy_array2) 
    plt.legend(["Local Floquet Hamiltonian","Random State"])
    plt.show()

def page_curve_for_random_conserved_case(N,NA,nks):
    entropy=0
    temp1=NA*np.log(2)
    temp2=0
    temp3=0
    temp4=0
    temp5=0
    for i in range(N):
        temp2 += -4*NA/N*nks[i]**2+8/3*NA/N*nks[i]**3-4/3*NA/N*nks[i]**4+NA/N*3/4
        temp3 += -4*NA**2/N**2*(nks[i]*(1-nks[i]))+8*NA**2/N**2*nks[i]*nks[i]*(1-nks[i])-16/3*NA**2/N**2*nks[i]**2*nks[i]*(1-nks[i])-8/3*NA**2/N**2*(nks[i]*(1-nks[i]))**2
        temp4 += -8/3*NA**3/N**3*(nks[i]*(1-nks[i]))**2
        #temp5 += -64/5*NA**4/N**4*(nks[i]*(1-nks[i]))**3+4/3*NA**4/N**4*(nks[i]*(1-nks[i]))**2
        temp5 += 4/3*NA**4/N**4*(nks[i]*(1-nks[i]))**2
    entropy=temp1+temp2+temp3+temp4+temp5
    quasi_particle_entropy=0
    for i in nks:
        if i>=1 or i<=0:
            continue
        quasi_particle_entropy+=-(1-i)*np.log2(1-i)-i*np.log2(i)
    quasi_particle_entropy=(NA/N-NA**2/N**2)*quasi_particle_entropy
    return entropy/np.log(2),quasi_particle_entropy

def subsystem_term_entropy_for_period2(N,NA):
    CM=generate_covariance_matrix_density_wave(2,N)
    U=transform_matrix_from_real_basis_to_conserved_basis(N)
    CM=U.dot(CM.dot(U.T.conj()))
    nks=np.diag(CM)
    nks=np.real(nks)
    return page_curve_for_random_conserved_case(N,NA,nks)

def calculate_power_two(NA,N):
    return NA**2/N

def calculate_power_four(NA,N):
    return 2*NA**3/N**2-NA**4/N**3

def calculate_power_six(NA,N):
    temp=(5+1/2)*NA**4/N**3-8*NA**5/N**4+4*NA**6/N**5
    temp2=5*NA**4/N**3-6*NA**5/N**4+2*NA**6/N**5
    '''
    temp=NA**6/N**5
    temp+=9/N**5*NA**5*(N-NA)
    if NA<=N/2:
        temp+=6/N**5*NA**4*(N-NA)*(N-2*NA)
    '''
    return temp,temp2

def calculate_power_ten(NA,N):
    #if NA>N/2:
    #    temp=NA**10/N**9+125*NA**9/N**9*(N-NA)
    #    return temp
    temp=120*NA**6/N**5-600*NA**7/N**6+1250*NA**8/N**7
    temp=temp-1225*NA**9/N**8+456*NA**10/N**9
    return temp

def calculate_power_eight(NA,N):
    temp=24*NA**5/N**4-72*NA**6/N**5+82*NA**7/N**6-33*NA**8/N**7
    return temp


def calculate_entropy_directly(NA,N):
    temp1=NA*np.log(2)
    temp2=-1/2*NA**2/N
    temp3=-1/6*NA**3/N**2
    temp4=-1/10*NA**4/N**3
    temp5=-0.06*NA**5/N**4
    temp6=-19/105*NA**6/N**5
    temp7=437/84*NA**7/N**6
    temp8=-6703/504*NA**8/N**7
    temp9=245/18*NA**9/N**8
    temp10=-76/15*NA**10/N**9
   # return (temp1+temp2+temp3+temp4+temp5+temp6+temp7+temp8+temp9+temp10)/np.log(2)
    return (temp1+temp2+temp3+temp4)/np.log(2)


def local_conservable_drawing():
    N_sites=input("N_sites:")
    N_sites=int(N_sites)
    N_size=N_sites
    Samplingtime=input("Sampling Time:")
    fileindex=input("FileIndex:")
    fs=open("./not_important_data/{}Nsite{}_Samplingtime{}.txt".format(fileindex,N_sites,Samplingtime),'r')
    data=fs.read().split()
    fs.close()
    Calculate_list=data[0::2]
    value=data[1::2]
    Calculate_list=[int(i) for i in Calculate_list]
    value=[float(i) for i in value]
    entropy_array2=[]
    entropy_array3=[]
    for j in Calculate_list:
        instance1=trun_k_probability(N_sites,N_sites//2,j)
        entropy_array2.append(instance1.direct_entropy_without_sum())
        j1=min(N_sites-j,j)
        entropy_array3.append(calculate_entropy_directly(j1, N_size))
    plt.rcParams.update({'font.size':30})
    plt.xlabel("Subsystem Size")
    plt.ylabel("Entropy")
    plt.plot(Calculate_list,value,Calculate_list,entropy_array2,Calculate_list,entropy_array3) 
    plt.legend(["Local Floquet Hamiltonian","Random State","Taylor Expansion"])
    plt.show()
def read_page_curve_comparison_continous():
    N_sites=input("N_sites:")
    N_sites=int(N_sites)
    Totaltime=input("Totaltime:")
    Samplingtime=input("Sampling Time:")
    T_period=input("T_period:")
    fs=open("./plot_data/Continous_Nsite{}_Tperiod{}Totaltime{}Samplingtime{}Create_byC.txt".format(N_sites,T_period,Totaltime,Samplingtime),"r")
    data=fs.read().split()
    fs.close()
    Calculate_list=data[0::2]
    value=data[1::2]
    Calculate_list=[int(i) for i in Calculate_list]
    value=[float(i) for i in value]
    entropy_array2=[]
    for j in Calculate_list:
        instance1=trun_k_probability(N_sites,N_sites//2,j)
        entropy_array2.append(instance1.direct_entropy_without_sum())
    plt.plot(Calculate_list,value,Calculate_list,entropy_array2) 
    plt.legend(["Local Floquet Hamiltonian","Random State"])
    plt.show()

def read_from_MF():
    N_sites=input("N_sites:")
    N_sites=int(N_sites)
    Totaltime=input("Totaltime:")
    Samplingtime=input("Sampling Time:")
    T_period=input("T_period:")
    method=input("Method:")
    fs=open("./plot_data/Nsite{}_Tperiod{}_Totaltime{}_Method{}MF.txt".format(N_sites,T_period,Totaltime,method),"r")
    data=fs.read().split()
    fs.close()
    Calculate_list=data[0::2]
    value=data[1::2]
    Calculate_list=[int(i) for i in Calculate_list]
    value=[float(i) for i in value]
    entropy_array2=[]
    for j in Calculate_list:
        instance1=trun_k_probability(N_sites,N_sites//2,j)
        entropy_array2.append(instance1.direct_entropy_without_sum())
    plt.rcParams.update({'font.size':30})
    plt.xlabel("Subsystem Size")
    plt.ylabel("Entropy")
    plt.plot(Calculate_list,value,Calculate_list,entropy_array2) 
    plt.legend(["Local Floquet Hamiltonian","Random State"])
    plt.show()


def Non_MF():
    N_sites=input("N_sites:")
    N_sites=int(N_sites)
    Totaltime=input("Totaltime:")
    Samplingtime=input("Sampling Time:")
    T_period=input("T_period:")
    method=input("Method:")
    fs=open("./plot_data/Nsite{}_Tperiod{}_Totaltime{}_Method{}Samplingtime{}.txt".format(N_sites,T_period,Totaltime,method,Samplingtime),"r")
    data=fs.read().split()
    fs.close()
    Calculate_list=data[0::2]
    value=data[1::2]
    Calculate_list=[int(i) for i in Calculate_list]
    value=[float(i) for i in value]
    entropy_array2=[]
    for j in Calculate_list:
        instance1=trun_k_probability(N_sites,N_sites//2,j)
        entropy_array2.append(instance1.direct_entropy_without_sum())
    plt.plot(Calculate_list,value,Calculate_list,entropy_array2) 
    plt.legend(["Local Floquet Hamiltonian","Random State"])
    plt.show()


def read_from_Hamiltonian_evolution():
    N_sites=input("N_sites:")
    N_sites=int(N_sites)
    N_size=N_sites
    Samplingtime=input("Sampling Time:")
    T_period=input("T_period:")
    Initialtime=input("Initial Time:")
    space=input("Space:")
    fileindex=input("Fileindex:")
    whether_random=input("Whether Random:")
    whether_random=int(whether_random)
    fs=open("./plot_data/Hamiltonian_evolution_Nsite{}_Tperiod{}_Initialtime{}_Samplingtime{}Spacing{}{}.txt".format(N_sites,T_period,Initialtime,Samplingtime,space,fileindex),"r")
    #fs=open("./plot_data/Nsize400norder0Samplingsize1000000Initialtime100000000.000000timestep0.200000.txt","r")
    data=fs.read().split()
    fs.close()
    Calculate_list=data[0::2]
    value=data[1::2]
    Calculate_list=[int(i) for i in Calculate_list]
    value=[float(i) for i in value]
    entropy_array2=[]
    entropy_array3=[]
    entropy_array4=[]
    for j in Calculate_list:
        instance1=trun_k_probability(N_sites,N_sites//int(space),j)
        entropy_array2.append(instance1.direct_entropy_without_sum())
        j1=min(j,N_sites-j)
        if whether_random==0:
            entropy_array3.append(calculate_entropy_directly(j1, N_size))
        if whether_random==1:
            entropy1,entropy2=subsystem_term_entropy_for_period2(N_size, j1)
            entropy_array3.append(entropy1)
            entropy_array4.append(entropy2)
    plt.rcParams.update({'font.size':40})
   # plt.xlabel(r"$f=N_A/N$",fontsize=40)
   # plt.ylabel(r"$\frac{\overline{S_A}}{N}$",fontsize=40)
   # plt.tight_layout()
    Calculate_list=np.array(Calculate_list)/N_size
    if whether_random==1:
#        plt.plot(Calculate_list,value,Calculate_list,entropy_array2,Calculate_list,entropy_array3,Calculate_list,entropy_array4) 
#       plt.legend(["Dynamical Page Curve","Random Gaussian State Page Curve","Taylor Expansion Result","Quasi Particle Picture"])
        plt.figure(figsize=[12,8],dpi=300)
        plt.plot(Calculate_list,np.array(value)/N_size,'b',marker='s')
        plt.plot(Calculate_list,np.array(entropy_array2)/N_size,'r',marker='d')
        plt.plot(Calculate_list,np.array(entropy_array3)/N_size,'g',marker='o')
        plt.tight_layout()
        plt.xlim(0.231,0.239)
        plt.xticks([0.232,0.234,0.236,0.238])
        plt.ylim(0.15,0.195)
        plt.yticks([0.15,0.16,0.17,0.18,0.19])
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)
        plt.savefig("../script/figs/{}Taylor_detailed.png".format(fileindex))
        plt.show()
        plt.figure(figsize=[12,8],dpi=300)
        plt.plot(Calculate_list,np.array(value)/N_size,'b',marker='s')
        plt.plot(Calculate_list,np.array(entropy_array2)/N_size,'r',marker='d')
        plt.plot(Calculate_list,np.array(entropy_array3)/N_size,'g',marker='o')
        plt.xlabel(r"$f=N_A/N$",fontsize=40)
        plt.ylabel(r"$\frac{\overline{S_A}}{N}$",fontsize=40)
        plt.xticks([0.0,0.2,0.4,0.6,0.8,1.0])
        plt.yticks([0.00,0.05,0.10,0.15,0.20,0.25])
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)
        plt.tight_layout()
        plt.legend(["Dynamical","Fully Random","Theory"],loc='upper right',prop={'size':30})
        plt.savefig("../script/figs/{}.png".format(fileindex))
      #  plt.show()
    if whether_random==0:
        #plt.plot(Calculate_list,value,'b',Calculate_list,entropy_array2,'r',Calculate_list,entropy_array3,'g')
        plt.figure(figsize=[12,8],dpi=300)
        plt.plot(Calculate_list,np.array(value)/N_size,'b',marker='s')
        plt.plot(Calculate_list,np.array(entropy_array2)/N_size,'r',marker='d')
        plt.plot(Calculate_list,np.array(entropy_array3)/N_size,'g',marker='o')
        plt.xlim(0.352,0.360)
        plt.ylim(0.2485,0.253)
        plt.tight_layout()
        plt.xticks([0.352,0.354,0.356,0.358,0.360])
        plt.yticks([0.249,0.250,0.251,0.252,0.253])
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)
        plt.savefig("../script/figs/{}Taylor_detailed.png".format(fileindex))
        plt.show()
        plt.figure(figsize=[12,8],dpi=300)
        plt.plot(Calculate_list,np.array(value)/N_size,'b',marker='s')
        plt.plot(Calculate_list,np.array(entropy_array2)/N_size,'r',marker='d')
        plt.plot(Calculate_list,np.array(entropy_array3)/N_size,'g',marker='o')
        plt.xlabel(r"$f=N_A/N$",fontsize=40)
        plt.ylabel(r"$\frac{\overline{S_A}}{N}$",fontsize=40)
        plt.xticks([0.0,0.2,0.4,0.6,0.8,1.0])
        plt.yticks([0.00,0.05,0.10,0.15,0.20,0.25])
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)
        plt.tight_layout()
        plt.legend(["Dynamical","Fully Random","Theory"],loc='upper right',prop={'size':30})
        plt.savefig("../script/figs/{}Taylor.png".format(fileindex))
        plt.show()

def draw_power_covariance_matrix():
    Nsize=input("Nsize:")
    norder=input("norder:")
    SamplingSize=input("SamplingSize:")    
    Initialtime=input("Initialtime:")
    timestep=input("timestep:")
    fs=open("./plot_data/Nsize{}norder{}Samplingsize{}Initialtime{}timestep{}.txt".format(Nsize,norder,SamplingSize,Initialtime,timestep),'r')
    data=fs.read().split()
    values_dynamical=[float(i) for i in data[1::2]]
    NASize=[float(i) for i in data[0::2]]
    values_theory=[]
    values_theory2=[]
    Nsize=int(Nsize)
    for j in range(Nsize):
        values_dynamical[j]/=Nsize
        if norder=='5':
            temp=calculate_power_ten(NASize[j],Nsize)
        elif norder=='4':
            temp=calculate_power_eight(NASize[j],Nsize)
        elif norder=='3':
            temp,temp2=calculate_power_six(NASize[j],Nsize)
            values_theory2.append(temp2/Nsize)
        elif norder=='2':
            temp=calculate_power_four(NASize[j],Nsize)
        elif norder=='1':
            temp=calculate_power_two(NASize[j],Nsize)
        else:
            print("order not defined")
            exit(0)
        values_theory.append(temp/Nsize)
        NASize[j]/=Nsize
    #plt.figure(figsize=[16.18,10],dpi=300)
    plt.plot(NASize,values_dynamical,'b',marker='s')
    plt.plot(NASize,values_theory,'g',marker='o')
    plt.xlabel(r"$f=N_A/N$",fontsize=40)
    if norder=='5':
        plt.ylabel(r"$\frac{\mathrm{Tr}\overline{X_A^10}}{N}$",fontsize=40)
    elif norder=='4':
        plt.ylabel(r"$\frac{\mathrm{Tr}\overline{X_A^8}}{N}$",fontsize=40)
    elif norder=='3':
        plt.plot(NASize,values_theory2,'r',marker='d')
        plt.ylabel(r"$\frac{\mathrm{Tr}\overline{X_A^6}}{N}$",fontsize=40)
    elif norder=='2':
        plt.ylabel(r"$\frac{\mathrm{Tr}\overline{X_A^4}}{N}$",fontsize=40)
    elif norder=='1':
        plt.ylabel(r"$\frac{\mathrm{Tr}\overline{X_A^2}}{N}$",fontsize=40)
    else:
        print("not defined order")
        exit(0)
    plt.xticks([0.0,0.2,0.4,0.6,0.8,1.0])
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.tight_layout()
    plt.legend(["Dynamical","Theory"],prop={'size':30})
    plt.show()
    


if __name__=="__main__":
    method=input("Please input method, 1 for Mf, 2 for continous, 3 for python_MF, 4 for non-MF, 5 for Hamiltonian evolution, 6 for local operator conserved case, 7 for power of covariance matrix:")
    if(method=='1'):
        read_page_curve_comparison()
    if(method=='2'):
        read_page_curve_comparison_continous()
    if(method=='3'):
        read_from_MF()
    if(method=='4'):
        Non_MF()
    if(method=='5'):
        read_from_Hamiltonian_evolution()
    if(method=='6'):
        local_conservable_drawing()
    if(method=='7'):
        draw_power_covariance_matrix()
