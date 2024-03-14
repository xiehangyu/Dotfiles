import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import sys

# Custom colormap
colors = [
    (1, 1, 1),  # White
    (0.9, 0.9, 0.9),  # Light gray
    (0.4, 0.4, 0.4),  # Dark gray
    (0, 0, 0.5),  # Dark blue
    (0, 0, 1),  # Blue
]
colors2=[
    (1, 1, 1),  # White
    (0.9, 0.9, 0.9),  # Light gray
    (0.4, 0.4, 0.4),  # Dark gray
    (0.5, 0, 0),  # Dark red
    (1, 0, 0),  # Red
]
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
cmap2=LinearSegmentedColormap.from_list("custom_cmap",colors2)
# Custom normalization function
class CustomNormalize(Normalize):
    def __init__(self,threshold=0.13,vmin=None, vmax=None, clip=None, normvalue=1.0):
        self.threshold = threshold
        self.norm=normvalue
        super().__init__(vmin, vmax, clip)
    def __call__(self, value, clip=None):
        x, y = [0, self.threshold, 1], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value/self.norm, x, y))




# Plotting
def plot_average_SFF():
    data1=open("./ensemble_average_twoCNOT_3_1round.txt","r")
    data1=data1.read().split()
    x=np.array([float(i) for i in data1[0::2]])
    y=np.array([float(i) for i in data1[1::2]])
    y=y-x
    plt.xlabel("Time t",fontsize=30)
    plt.ylabel(r"$\langle SFF\rangle -t$",fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.plot(x,y)
    data2=open("./ensemble_average_twoCNOT_3_2round.txt","r")
    data2=data2.read().split()
    x=np.array([float(i) for i in data2[0::2]])
    y=np.array([float(i) for i in data2[1::2]])
    y=y-x
    plt.plot(x,y)
    data3=open("./ensemble_average_twoCNOT_4_1round.txt","r")
    data3=data3.read().split()
    x=np.array([float(i) for i in data3[0::2]])
    y=np.array([float(i) for i in data3[1::2]])
    y=y-x
    plt.plot(x,y)
    data4=open("./ensemble_average_twoCNOT_4_2round.txt","r")
    data4=data4.read().split()
    x=np.array([float(i) for i in data4[0::2]])
    y=np.array([float(i) for i in data4[1::2]])
    y=y-x
    plt.plot(x,y)
    data5=open("./ensemble_average_twoCNOT_5_1round.txt","r")
    data5=data5.read().split()
    x=np.array([float(i) for i in data5[0::2]])
    y=np.array([float(i) for i in data5[1::2]])
    y=y-x
    plt.plot(x,y)
    data6=open("./ensemble_average_twoCNOT_5_2round.txt","r")
    data6=data6.read().split()
    x=np.array([float(i) for i in data6[0::2]])
    y=np.array([float(i) for i in data6[1::2]])
    y=y-x
    plt.plot(x,y)
    data7=open("./ensemble_average_twoCNOT_6_1round.txt","r")
    data7=data7.read().split()
    x=np.array([float(i) for i in data7[0::2]])
    y=np.array([float(i) for i in data7[1::2]])
    y=y-x
    plt.plot(x,y)
    data8=open("./ensemble_average_twoCNOT_6_2round.txt","r")
    data8=data8.read().split()
    x=np.array([float(i) for i in data8[0::2]])
    y=np.array([float(i) for i in data8[1::2]])
    y=y-x
    plt.plot(x,y)
    plt.legend([r"$N=3,1$",r"$N=3,2$",r"$N=4,1$",r"$N=4,2$",r"$N=5,1$",r"$N=5,2$",r"$N=6,1$",r"$N=6,2$"],fontsize=40)
    plt.show()


def plot_heatmap_correlation2(filename,colors='red',thresheld=0.13):

# Plotting
    if colors=='red':
        usemap=cmap2
    if colors=='blue':
        usemap=cmap
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams["font.family"] = "STIXGeneral"
    plt.figsize=(10,5)
    plt.rcParams.update({'font.size': 40})
    fp=open(filename,"r") 
    data=fp.read().split()
    t=data[0::3]
    x=data[1::3]
    z=data[2::3]
    tsize=len(set(t))
    xsize=len(set(x))
    t=[float(i) for i in t]
    x=[float(i) for i in x]
    z=[abs(complex(i)) for i in z]
    norm2=max(z)
    z=np.reshape(z,(tsize,xsize))
    norm=CustomNormalize(threshold=thresheld,normvalue=norm2)
    # Set number of squares in x and y direction
    
    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(12,8))
    
    # Draw the heatmap
    heatmap = ax.imshow(z,extent=[min(x),max(x),max(t),min(t)], cmap=usemap,norm=norm,origin='upper')
    print(z)
    
    
    # Label x and y ticks with the corresponding row/column labels
    #ax.set_xticklabels(np.arange(0, 14))
    #ax.set_yticklabels(np.arange(0, 6))
    
    # Set axis labels
    ax.set_xlabel(r'Location $2i$',fontsize=40)
    ax.set_ylabel(r'Time $t$',fontsize=40)
    ax.tick_params(axis='x', labelsize=40)
    ax.tick_params(axis='y', labelsize=40)
    ax.set_xticks([0,10,20,30,40])
    ax.set_yticks([0,5,10,15,20])
    ax.tick_params(axis='both',direction='in')
    ax.set_ylim([1,np.max(ax.get_ylim())])
    # Add a colorbar
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=usemap), ax=ax,shrink=0.5)
    cbar.set_label(r"$|C_{ij}(t)|$",fontsize=40)
    cbar.set_ticks([0,0.5,1.0])
    plt.tight_layout()
    plt.savefig("/home/xiehangyu/yxh/study/mpq/quantum_circuit/script/Figs/"+filename.replace(".txt",".pdf"),dpi=300)
    # Show the plot
    plt.show()

def plot_average_SFF():
    data1=open("./ensemble_average_twoCNOT_3_1round.txt","r")
    data1=data1.read().split()
    x=np.array([float(i) for i in data1[0::2]])
    y=np.array([float(i) for i in data1[1::2]])
    y=y-x
    plt.xlabel("Time t",fontsize=30)
    plt.ylabel(r"$\langle SFF\rangle -t$",fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.plot(x,y)
    data2=open("./ensemble_average_twoCNOT_3_2round.txt","r")
    data2=data2.read().split()
    x=np.array([float(i) for i in data2[0::2]])
    y=np.array([float(i) for i in data2[1::2]])
    y=y-x
    plt.plot(x,y)
    data3=open("./ensemble_average_twoCNOT_4_1round.txt","r")
    data3=data3.read().split()
    x=np.array([float(i) for i in data3[0::2]])
    y=np.array([float(i) for i in data3[1::2]])
    y=y-x
    plt.plot(x,y)
    data4=open("./ensemble_average_twoCNOT_4_2round.txt","r")
    data4=data4.read().split()
    x=np.array([float(i) for i in data4[0::2]])
    y=np.array([float(i) for i in data4[1::2]])
    y=y-x
    plt.plot(x,y)
    data5=open("./ensemble_average_twoCNOT_5_1round.txt","r")
    data5=data5.read().split()
    x=np.array([float(i) for i in data5[0::2]])
    y=np.array([float(i) for i in data5[1::2]])
    y=y-x
    plt.plot(x,y)
    data6=open("./ensemble_average_twoCNOT_5_2round.txt","r")
    data6=data6.read().split()
    x=np.array([float(i) for i in data6[0::2]])
    y=np.array([float(i) for i in data6[1::2]])
    y=y-x
    plt.plot(x,y)
    data7=open("./ensemble_average_twoCNOT_6_1round.txt","r")
    data7=data7.read().split()
    x=np.array([float(i) for i in data7[0::2]])
    y=np.array([float(i) for i in data7[1::2]])
    y=y-x
    plt.plot(x,y)
    data8=open("./ensemble_average_twoCNOT_6_2round.txt","r")
    data8=data8.read().split()
    x=np.array([float(i) for i in data8[0::2]])
    y=np.array([float(i) for i in data8[1::2]])
    y=y-x
    plt.plot(x,y)
    plt.legend([r"$N=3,1$",r"$N=3,2$",r"$N=4,1$",r"$N=4,2$",r"$N=5,1$",r"$N=5,2$",r"$N=6,1$",r"$N=6,2$"],fontsize=30)
    plt.show()

def plot_heatmap_correlation():

# Plotting

    fp=open("./correlation_function_randomU2.txt","r") 
    data=fp.read().split()
    x=data[4::3]
    y=data[5::3]
    z=data[6::3]
    x=[float(i) for i in x]
    y=[float(i) for i in y]
    z=[abs(complex(i)) for i in z]
    z=np.reshape(z,(6,14))
    # Set number of squares in x and y direction
    
    # Create a figure and axis object
    fig, ax = plt.subplots()
    
    # Draw the heatmap
    heatmap = ax.imshow(z/1000.0, cmap=cmap,norm=norm)
    
    
    # Label x and y ticks with the corresponding row/column labels
    #ax.set_xticklabels(np.arange(0, 14))
    #ax.set_yticklabels(np.arange(0, 6))
    
    # Set axis labels
    ax.set_xlabel('Site Number',fontsize=30)
    ax.set_ylabel('Time t',fontsize=30)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_ylim(ax.get_ylim()[::-1])
    # Add a colorbar
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    
    # Show the plot
    plt.show()

def linear_function1(x,a,b):
    return a*x+b

def fitting_quench():
    data=open("./thermalize_quantum_quench.txt","r")
    data2=open("./thermalize_quantum_quenchrandom2.txt",'r')
    data3=open("./thermalize_quantum_quenchanother.txt","r")
    data=data.read().split()
    data2=data2.read().split()
    data3=data3.read().split()
    plt.figure(figsize=(15,10))
    plt.rcParams.update({'font.size': 40})
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams["font.family"] = "STIXGeneral"
    x=np.array([float(i) for i in data[0::4]])/2
    y=np.array([abs(complex(i)) for i in data[1::4]])
    y=np.log10(y)
    x2=np.array([float(i) for i in data2[0::4]])/2
    y2=np.array([abs(complex(i)) for i in data2[1::4]])
    y2=np.log10(y2)
    x3=np.array([float(i) for i in data3[0::4]])/2
    y3=np.array([abs(complex(i)) for i in data3[1::4]])
    plt.scatter(x,y,s=60,c='b')
    plt.scatter(x2,y2,s=60,c='g')
    plt.scatter(x3,y3,s=60,c='r')
    plt.legend(["Random1","Random2","Non Ergodic","","",""],fontsize=40)
    xs=np.linspace(min(x),max(x),1000)
    popt,pcov=curve_fit(linear_function1,x,y)
    xs2=np.linspace(min(x2),max(x2),1000)
    popt2,pcov2=curve_fit(linear_function1,x2,y2)
    xs3=np.linspace(min(x3),max(x3),1000)
    popt3,pcov3=curve_fit(linear_function1,x3,y3)
    plt.plot(xs,linear_function1(xs,*popt),'b',linewidth=3)
    plt.plot(xs2,linear_function1(xs2,*popt2),'g',linewidth=3)
    plt.plot(xs3,linear_function1(xs3,*popt3),'r',linewidth=3)
    plt.tick_params(axis='both',direction='in')
    plt.xlabel("Time t",fontsize=40)
    plt.ylabel(r"$\log_{10}|\langle\psi|\hat{O}(t)|\psi\rangle|$",fontsize=40)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.tight_layout()
    plt.savefig("/home/xiehangyu/yxh/study/mpq/quantum_circuit/script/Figs/thermal_quantum_quench.pdf",dpi=300)
    plt.show()

def fittingplot():
    data=open("./decay_of_correlationfunction.txt","r")
    data=data.read().split()
    plt.figure(figsize=(12,8))
    plt.rcParams.update({'font.size': 40})
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams["font.family"] = "STIXGeneral"
    x=np.array([float(i) for i in data[0::3]])/2
    y1=np.array([float(i) for i in data[1::3]])
    y2=np.array([float(i) for i in data[2::3]])
    x=x[0:60]
    y1=y1[0:60]
    y2=y2[0:60]
    y1=np.log10(y1)
    y2=np.log10(y2)
    plt.scatter(x,y1,s=60,c='b')
    plt.scatter(x,y2,s=60,c='g')
    plt.legend([r"i=j",r"i=j+t","",""],fontsize=40)
    xs=np.linspace(min(x),max(x),1000)
    popt1,pcov1=curve_fit(linear_function1,x,y1)
    popt2,pcov2=curve_fit(linear_function1,x,y2)
    plt.plot(xs,linear_function1(xs,*popt1),'b',linewidth=3)
    plt.plot(xs,linear_function1(xs,*popt2),'g',linewidth=3)
    plt.tick_params(axis='both',direction='in')
    plt.xlabel("Time t",fontsize=40)
    plt.ylabel(r"$\log_{10}C_{ij}(t)})$",fontsize=40)
    plt.xticks(fontsize=40)
    plt.yticks([0,-10,-20,-30])
    plt.yticks(fontsize=40)
    plt.tight_layout()
    plt.savefig("/home/xiehangyu/yxh/study/mpq/quantum_circuit/script/Figs/decay_of_correlationfunction.pdf",dpi=300)
    plt.show()
    
def plot_heatmap_correlation3(filename):
    data=open(filename,"r")
    data=data.read().split()
    plt.rcParams.update({'font.size': 40})
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams["font.family"] = "STIXGeneral"
    t=np.array([float(i) for i in data[0::3]])
    x=np.array([float(i) for i in data[1::3]])
    z=np.array([abs(complex(i)) for i in data[2::3]])
    norm2=np.max(z)
    z[z<=1E-20]=0
    orig=np.copy(z)
    z=np.log10(z)
    z[orig==0]=np.nan
    tsize=len(set(t))
    xsize=len(set(x))
    z=np.reshape(z,(tsize,xsize))
    # Set number of squares in x and y direction
    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(12,8))
    # Draw the heatmap
    cmap=plt.cm.viridis.copy()
    cmap.set_bad(color='white')
    heatmap = ax.imshow(z,extent=[min(x),max(x),max(t),min(t)], cmap=cmap,vmax=0,vmin=-15)
    # Label x and y ticks with the corresponding row/column labels
    #ax.set_xticklabels(np.arange(0, 14))
    #ax.set_yticklabels(np.arange(0, 6))
    # Set axis labels
    ax.set_xlabel(r'Location $2i$',fontsize=40)
    ax.set_ylabel(r'Time $t$',fontsize=40)
    ax.set_xticks([0,10,20,30,40])
    ax.set_yticks([0,5,10,15,20])
    ax.tick_params(axis='x', labelsize=40)
    ax.tick_params(axis='y', labelsize=40)
    ax.tick_params(axis='both',direction='in')
    ax.set_ylim([1,np.max(ax.get_ylim())])
    # Add a colorbar
    cbar = fig.colorbar(heatmap,ax=ax,shrink=0.5)
    cbar.set_label(r"$\log_{10}|C_{ij}(t)|$",fontsize=40)
    plt.tight_layout()
    plt.savefig("/home/xiehangyu/yxh/study/mpq/quantum_circuit/script/Figs/"+filename.replace(".txt",".pdf"),dpi=300)
    # Show the plot
    plt.show()



    


if __name__=='__main__':
    if sys.argv[1]=='1':
        plot_heatmap_correlation2("correlation_function_support2.txt",'blue',0.13)
    if sys.argv[1]=='2':
        plot_heatmap_correlation3("correlation_function_support3.txt")
    if sys.argv[1]=='3':
        fittingplot()
    if sys.argv[1]=='4':
        fitting_quench()