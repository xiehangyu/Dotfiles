import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
frame_width=2
def rising_function(x,I0,tau,d):
    return I0*(1-np.exp(-np.power(x,d)/tau))

def falling_function(x,A,tau1,B,tau2):
    return A*np.exp(-x/tau1)+B*np.exp(-x/tau2)
def FigureC():
    plt.figure(figsize=(10,8))
    fp1=open('./400nm_vaccum_excelent.txt','r')
    fp2=open('./400nm_air_excelent.txt','r')
    data1=fp1.read().split()
    data2=fp2.read().split()
    x1=data1[4::4]
    y1=data1[5::4]
    x2=data2[4::4]
    y2=data2[5::4]
    x1=[float(i) for i in x1]
    y1=[float(i)*1E11 for i in y1]
    x2=[float(i) for i in x2]
    y2=[float(i)*1E11 for i in y2]
    xdraw1=x1[1083:1117]
    xdraw2=x2[1081:1118]
    ydraw1=y1[1083:1117]
    ydraw2=y2[1081:1118]
    ydraw2=np.array(ydraw2)
    xdraw1=np.array(xdraw1)
    xdraw2=np.array(xdraw2)
    ydraw1=np.array(ydraw1)
    xdraw2=xdraw2-xdraw2[0]
    xdraw1=xdraw1-xdraw1[0]
    print(xdraw2)
    ydraw1=ydraw1-6.657
    ydraw2=ydraw2-10.46
    print(ydraw2)
    print(falling_function(xdraw2[2],5.428,0.8324,16.84,0.05752))
    xmin1=min(xdraw1)
    xmax1=max(xdraw1)
    xmin2=min(xdraw2)
    xmax2=max(xdraw2)
    xs1=np.linspace(xmin1,xmax1,500)
    xs2=np.linspace(xmin2,xmax2,500)
    yfit1=[]
    for i in xs1:
        yfit1.append(falling_function(i,13.16,0.0575,5.449,0.698))
    yfit2=[]
    for i in xs2:
        yfit2.append(falling_function(i,5.428,0.8324,16.84,0.05752))
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.plot(xs1,yfit1,'r',linewidth=3)
    plt.plot(xs2,yfit2,'b',linewidth=3)
    plt.scatter(xdraw1,ydraw1,s=150,marker='o', facecolors='none', edgecolors='red')
    plt.scatter(xdraw2,ydraw2,s=150,marker='o', facecolors='none', edgecolors='blue')
    plt.xlabel(r"Time $(s)$",fontsize=40)
    plt.ylabel(r"Current $(\times 10^{-2}nA)$",fontsize=40)
    plt.xticks(fontsize=40)
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.yticks(fontsize=40)
    plt.ylim([0,25])
    ax=plt.gca()
    ax.spines['bottom'].set_linewidth(frame_width)
    ax.spines['top'].set_linewidth(frame_width)
    ax.spines['left'].set_linewidth(frame_width)
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams['font.family'] = 'STIXGeneral'
    ax.spines['right'].set_linewidth(frame_width)
    plt.legend(["Fit Vacuum","Fit Air","Experiment Vacuum","Experiment Air"],fontsize=30)
    plt.tight_layout()
    plt.savefig("FigureC.pdf",dpi=300)
    plt.show()

def FigureD():
    plt.figure(figsize=(10,8))
    fp1=open('./400nm_vaccum_excelent.txt','r')
    fp2=open('./400nm_air_excelent.txt','r')
    data1=fp1.read().split()
    data2=fp2.read().split()
    x1=data1[4::4]
    y1=data1[5::4]
    x2=data2[4::4]
    y2=data2[5::4]
    x1=[float(i) for i in x1]
    y1=[float(i)*1E11 for i in y1]
    x2=[float(i) for i in x2]
    y2=[float(i)*1E11 for i in y2]
    xdraw1=x1[605:642]
    xdraw2=x2[601:638]
    ydraw1=y1[605:642]
    ydraw2=y2[601:638]
    ydraw2=np.array(ydraw2)
    xdraw1=np.array(xdraw1)
    xdraw2=np.array(xdraw2)
    ydraw1=np.array(ydraw1)
    xdraw2=xdraw2-xdraw2[0]
    xdraw1=xdraw1-xdraw1[0]
    ydraw1=ydraw1-6.843
    ydraw2=ydraw2-10.38
    xmin1=min(xdraw1)
    xmax1=max(xdraw1)
    xmin2=min(xdraw2)
    xmax2=max(xdraw2)
    xs1=np.linspace(xmin1,xmax1,500)
    xs2=np.linspace(xmin2,xmax2,500)
    yfit1=[]
    for i in xs1:
        yfit1.append(rising_function(i,18.74,0.2545,0.5026))
    yfit2=[]
    for i in xs2:
        yfit2.append(rising_function(i,24.44,0.2686,0.4591))
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.plot(xs1,yfit1,'r',linewidth=3)
    plt.plot(xs2,yfit2,'b',linewidth=3)
    plt.ylim([0,25])
    plt.scatter(xdraw1,ydraw1,s=150,marker='o', facecolors='none', edgecolors='red')
    plt.scatter(xdraw2,ydraw2,s=150,marker='o', facecolors='none', edgecolors='blue')
    plt.xlabel(r"Time $(s)$",fontsize=40)
    plt.ylabel(r"Current $(\times 10^{-2}nA)$",fontsize=40)
    plt.xticks(fontsize=40)
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.yticks(fontsize=40)
    ax=plt.gca()
    ax.spines['bottom'].set_linewidth(frame_width)
    ax.spines['top'].set_linewidth(frame_width)
    ax.spines['left'].set_linewidth(frame_width)
    ax.spines['right'].set_linewidth(frame_width)
    plt.legend(["Fit Vacuum","Fit Air","Experiment Vacuum","Experiment Air"],fontsize=30)
    plt.tight_layout()
    plt.savefig("FigureD.pdf",dpi=300)
    plt.show()

def FigureB():
    plt.figure(figsize=(10,4))
    fp1=open('./400nm_vaccum_excelent.txt','r')
    fp2=open('./400nm_air_excelent.txt','r')
    data1=fp1.read().split()
    data2=fp2.read().split()
    x1=data1[4::4]
    x2=data2[4::4]
    y1=data1[5::4]
    y2=data2[5::4]
    x1=x1[455:]
    y1=y1[455:]
    x2=x2[455:]
    y2=y2[455:]
    x1=[float(i) for i in x1]
    y1=[float(i) for i in y1]
    x2=[float(i) for i in x2]
    y2=[float(i) for i in y2]
    x1=np.array(x1)
    y1=np.array(y1)
    x2=np.array(x2)
    y2=np.array(y2)
    temp=x1[0]
    x1=x1-x1[0]
    x2=x2-x2[0]
    y1=y1*1E11
    y2=y2*1E11
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams['font.family'] = 'STIXGeneral'
    ax=plt.gca()
    ax.spines['bottom'].set_linewidth(frame_width)
    ax.spines['top'].set_linewidth(frame_width)
    ax.spines['left'].set_linewidth(frame_width)
    ax.spines['right'].set_linewidth(frame_width)
    plt.plot(x1,y1,'r',linewidth=3)
    plt.plot(x2,y2,'b',linewidth=3)
    plt.xlabel(r"Time $(s)$",fontsize=30)
    plt.ylabel(r"Current $(\times 10^{-2}nA)$",fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.text(196-temp,30, 'Air', ha='center', va='bottom',fontsize=30)
    plt.text(196-temp,10, 'Vacuum', ha='center', va='top',fontsize=30)
    plt.tight_layout()
    plt.savefig("FigureB.pdf",dpi=300)
    plt.show()




if __name__=='__main__':
    FigureC()
    