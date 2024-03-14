from first_trial import *
import matplotlib.pyplot as plt

if __name__=='__main__':
    instance1=trun_k_probability(200,100,80)
    mutual_infos=[]
    subsys=[]
    for i in range(0,81):
        mutual_infos.append(instance1.average_mutual_information(i))
        subsys.append(i)
    plt.figure(figsize=[12,8],dpi=300)
    plt.plot(subsys,mutual_infos)
    plt.tight_layout()
    plt.xlabel("Subsystem Size",fontsize=30)
    plt.ylabel("Mutual Information",fontsize=30)
    plt.tick_params(labelsize=30)
    plt.show()