# -*- coding:gbk -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
frame_width = 2


def ReadData(filename):
    # filename = ['415_20_1.txt', '415_40_1.txt', '415_60_1.txt', '415_80_1.txt', '415_100.txt', '415_120.txt', '415_140.txt', '415_160.txt', '415_180.txt', '415_200.txt']
    num = len(filename)
    x = []
    y = []
    for i in range(num):
        fp = open(filename[i], 'r')
        data = fp.read().split()
        x_t = data[6::4]
        y_t = data[5::4]
        x_t = [float(i) for i in x_t]
        y_t = [float(i)*1E9 for i in y_t]
        x_t = np.array(x_t)
        y_t = np.array(y_t)
        x.append(x_t)
        y.append(y_t)
        '''
        建议为了安全起见，还是加上fp.close()吧。
        python和C不一样，如果文件不关闭的话，那么原先写入文件的内容就不会更新到硬盘上。这样的话，如果程序意外终端，
        就会导致原先写入的内容也丢失了。我之前就犯过这样的问题。
        '''
        fp.close()
    return x, y

def HT_Plot(x,y,colormap='hsv',figname="HT_AVLog1.pdf",T=None):
    plt.figure(figsize=(10,8))
    '''
    换了一个函数取colormap
    colors = mpl.cm.get_cmap(colormap, 15)
    c = np.linspace(0, 1, 10)
    c = colors(c)
    '''
    total_color_range=12
    c = plt.cm.Set1(np.linspace(0, 1, total_color_range))
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams['font.family'] = 'Arial'
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(frame_width)
    ax.spines['top'].set_linewidth(frame_width)
    ax.spines['left'].set_linewidth(frame_width)
    ax.spines['right'].set_linewidth(frame_width)
    num = len(x)
    for i in range(num):
        plt.plot(x[i], y[i], c=c[i], linewidth=2)
    '''
    在log作图时，建议还是加上minor ticks（或称subticks）吧。它的意思是，除了显示例如1E-1,1E0,1E1的坐标外，还会在指定的子坐标上（例如2，3，4，）上显示刻线。我感觉这点还是挺有用的，可以帮助在log绘图上更轻易地读数。
    '''
    plt.minorticks_on()
    '''
    设置minor ticks(或称subticks)的位置。
    '''
    subticks=np.array([i*0.1 for i in range(1,10)])
    subticks=list(subticks)+list(subticks*10)+list(subticks*100)
    '''
    在yscale中的symlog其实并不一直是按照对数坐标画图的。在靠近原点的附近(-linthresh,linthresh)是按照线性坐标画图的。linscale参数就指定了这段线性区间(0,linthresh)在图上的长度相对于对数坐标单个区间（例如区间10-100）的比例。如果linthresh与linscale设置的值不匹配， 就会使得坐标区间看起来不是均匀分布的。在python的default value中，这两者的设置是不匹配的。为了使得坐标区间看起来均匀分布，需要手动设置两者的值
    '''
    plt.yscale('symlog',linthresh=0.1,linscale=1.0,subs=subticks)
    '''
   在一般的科研文献作图中，我建议fontsize至少取到30吧，我一般取到40（做PPT时也一样），不然可能会在文献里看不清楚。
    '''
    plt.xlabel(r"Voltage $(V)$", fontsize=20,ha='right')
    plt.ylabel(r"Current $(nA)$", fontsize=20,labelpad=-20)
    plt.xticks(fontsize=20)
    plt.xticks([-10,-5,0,5,10])
    plt.yticks(fontsize=20)
    plt.xlim([-10, 10])
    plt.ylim([-30, 30])
    plt.tick_params(axis='x', direction='in',width=2)
    plt.tick_params(axis='y', direction='in',width=2)
    '''
    我建议还是legend吧，这样可以很方便得看出每一条曲线对应什么温度。
    '''
    plt.legend([r"$"+str(i)+"K$" for i in T], fontsize=20)
    '''
    使minor ticks的指向也为'in'
    '''
    plt.tick_params(axis='y', which='minor', direction='in',width=2)
    plt.tick_params(axis='x',which='minor',direction='in',width=2)
    # plt.text(205,30, 'Air', ha='center', va='bottom',fontsize=30)
    # plt.text(200,10, 'Vacuum', ha='center', va='top',fontsize=30)
    plt.tight_layout()
    plt.savefig(figname, dpi=300)
    plt.show()


if __name__ == '__main__':
    filename = ['415_20_1.txt', '415_40_1.txt', '415_60_1.txt', '415_80_1.txt', '415_100.txt', '415_120.txt',
                '415_140.txt', '415_160.txt', '415_180.txt', '415_200.txt']
    T = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0]
    T = [(i+273.15) for i in T]
    x, y = ReadData(filename)

    '''
    这儿是不是笔误了？我感觉这里传进去的参数应该是x,而不是T?
    HT_Plot(T, y, colormap='hsv',figname="HT_AVLog1.pdf")
    '''
    HT_Plot(x, y, colormap='hsv',figname="HT_AVLog1.pdf",T=T)
    '''
    加了一个T参数，用来写legend。
    '''
