import matplotlib.pyplot as plt
from scipy.special import digamma
import numpy as np
from first_trial import *
def digamma_function(n):
    return digamma(n)


def entropy_calculation(N,NA):
    if NA>N/2:
        NA=N-NA
    temp_value=digamma_function(2**N+1)
    temp_value=temp_value-digamma_function(2**(N-NA)+1)
    temp_value=temp_value-(2**(NA)-1)/(2**(N-NA+1))
    return temp_value/np.log(2)

def thermal_dynamical_entropy_density(f):
    if f>1/2:
        f=1-f
    return f

def draw_page_curve():
    N=10
    fs=[]
    entropy1_densitys=[]
    entropy2_densitys=[]
    for i in range(0,101):
        f=i*0.01
        fs.append(f)
        entropy1=entropy_calculation(N, f*N)
        entropy1_densitys.append(entropy1/N)
        entropy2_densitys.append(thermal_dynamical_entropy_density(f))
    plt.rcParams.update({"font.size":30})
    plt.plot(fs,entropy1_densitys,fs,entropy2_densitys)
    plt.legend([r"Page Curve for $N=10$",r"Page Curve for $N\to\infty$"])
    plt.xlabel(r"$f=N_A/N$")
    plt.ylabel(r"$\frac{<S_A>}{N}$")
    plt.show()

def entropy_density_page_curve(N):
    sites=[]
    entropys=[]
    for i in range(1,N):
        sites.append(i) 
        entropys.append(entropy_calculation(N, i)/N)
    plt.plot(sites,entropys) 
    plt.show()

def mutual_information(N,NA,NA1):
    NA2=NA-NA1
    entropy1=entropy_calculation(N, NA)
    entropy2=entropy_calculation(N, NA1)
    entropy3=entropy_calculation(N, NA2)
    return entropy2+entropy3-entropy1


def mutual_information_page(N,NA):
    sites=[0]
    mutual_infs=[0]
    for i in range(1,NA):
        sites.append(i)
        mutual_infs.append(mutual_information(N, NA, i))
    mutual_infs.append(0)
    sites.append(NA)
    plt.figure(figsize=[12,8],dpi=300)
    plt.plot(sites,mutual_infs)
    plt.xticks(fontsize=40)
    plt.xticks([0,20,40,60,80])
    plt.yticks(fontsize=40)
    plt.xlabel(r"$N_{A_{1}}$",fontsize=40)
    plt.ylabel(r"$\langle I(A_1:A_2)\rangle$",fontsize=40)
    plt.tight_layout()
    plt.savefig("../script/figs/mutual_information_N{}_NA{}_interacting.png".format(N,NA))
    plt.show()


def different_mutual_information_with_N(NA):
    sites=[]
    mutus=[]
    for N in range(2*NA,63):
        sites.append(N)
        mutus.append(mutual_information(N, NA, NA//2))
    fs=open("./plot_data/interaction_mutual_information_with_N.txt","w")
    for i in range(len(sites)):
        fs.write("{}\t{}\n".format(sites[i],mutus[i]))
    fs.close()


def drawing_for_CRFGE():
    N=20
    instance1=trun_k_probability(N, 10, 1)
    instance2=trun_k_probability(N, 5, 1)
    entropy1s=[]
    entropy2s=[]
    NAs=[]
    entropy3s=[]
    entropy4s=[]
    ms=[]
    instance3=trun_k_probability(N, 1, 10)
    instance4=trun_k_probability(N, 1, 5)
    for i in range(1,N):
        instance1.change_parameter(N, 10, i)
        instance2.change_parameter(N, 5, i)
        instance3.change_parameter(N, i, 10)
        instance4.change_parameter(N, i, 5)
        NAs.append(i)
        entropy1s.append(instance1.direct_entropy_without_sum()/N)
        entropy2s.append(instance2.direct_entropy_without_sum()/N)
        ms.append(i)
        entropy3s.append(instance3.direct_entropy_without_sum()/N)
        entropy4s.append(instance4.direct_entropy_without_sum()/N)
    plt.rcParams.update({"font.size":30}) 
    plt.plot(NAs,entropy1s,NAs,entropy2s)
    plt.xlabel("Subsystem Size")
    plt.ylabel(r"$\frac{<S_A>}{N}$")
    plt.legend(["Charge Number=10","Charge Number=5"])
    plt.show()
    plt.plot(ms,entropy3s,ms,entropy4s)
    plt.xlabel("Charge Number")
    plt.ylabel(r"$\frac{<S_A>}{N}$")
    plt.legend(["Subsystem Size=10","Subsystem Size=5"])
    plt.show()
    
def Gaussian_entropy(N,NA):
    temp_value=(N-0.5)*digamma(2*N)
    temp_value+=(0.5+NA-N)*digamma(2*N-2*NA)
    temp_value+=(0.25-NA)*digamma(N)
    temp_value-=0.25*digamma(N-NA)
    temp_value-=NA
    return temp_value/np.log(2)

def Gaussian_entropy_thermal_dynamical_limit(f):
    if f>0.5:
        f=1-f
    temp_value=f*(np.log(2)-1)
    temp_value+=(f-1)*np.log(1-f)
    return temp_value/np.log(2)

def draw_picture_comparison_gaussian_interaction():
    fs=[]
    entropys_in=[]
    entropys_ga=[]
    for i in np.linspace(0,0.5,50):
        fs.append(i)
        entropys_in.append(thermal_dynamical_entropy_density(i))
        entropys_ga.append(Gaussian_entropy_thermal_dynamical_limit(i))
    j=len(fs)-1
    while j>=0:
        fs.append(1-fs[j])
        entropys_in.append(entropys_in[j])
        entropys_ga.append(entropys_ga[j])
        j=j-1
    plt.figure(figsize=[16.18,10],dpi=300)
    plt.plot(fs,entropys_in,'k--')
    plt.plot(fs,entropys_ga,'b-')
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.annotate('Microscopic',xy=(0,0),fontsize=40,xytext=(0.08,0.2),ha='center')
    plt.annotate('',xy=(0,0),xytext=(0.076,0.19),arrowprops={'facecolor':'black'})
    plt.annotate('Macroscopic',xy=(0.38,0.3),fontsize=40,xytext=(0.3,0.45),ha='center')
    plt.annotate('',xy=(0.38,0.3),xytext=(0.304,0.442),arrowprops={'facecolor':'black'})
    plt.xlabel(r"$f=N_A/N$",fontsize=40)
    plt.ylabel(r"$\frac{\langle{S_A}\rangle}{N}$",fontsize=40)
    plt.legend(["Interacting Page Curve","Free-fermion Page Curve"],prop={'size':30})
    plt.tight_layout()
    plt.savefig("../script/figs/illustration_of_mi_or_ma.png")



if __name__=="__main__":
    #draw_picture_comparison_gaussian_interaction()
    mutual_information_page(50,20)