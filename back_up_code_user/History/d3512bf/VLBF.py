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
        为了程序安全起见，最好加上fp.close()吧。
        python和C不一样, python有一个毛病，就是如果文件不关闭，那么写进文件里的内容始终不会更新到磁盘上。尤其是在写文件的时候，如果程序意外中断，那么之前计划写进文件里的内容也会全部丢失。之前我就出过这样的问题。
        '''
    return x, y

def HT_Plot(x,y,colormap='hsv',figname="HT_AVLog1.pdf",T=None):
    plt.figure
    colors = mpl.cm.get_cmap(colormap, 10)
    c = np.linspace(0, 1, 10)
    c = colors(c)
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams['font.family'] = 'STIXGeneral'
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(frame_width)
    ax.spines['top'].set_linewidth(frame_width)
    ax.spines['left'].set_linewidth(frame_width)
    ax.spines['right'].set_linewidth(frame_width)
    num = len(x)
    for i in range(num):
        plt.plot(x[i], y[i], c=c[i], linewidth=2)
    '''
    在log坐标下，可以开启minor ticks,也就是除了1E-1,1E0,1E1外，还会将中间的坐标值，例如2，4，20，40等也显示出来。我感觉这在log坐标下还不是比较重要的，可以更方便读出数据的值。
    '''
    plt.minorticks_on()
    '''
    设置minorticks（或称做subticks）的数值，也就是在哪些点会标记出小坐标。
    '''
    subticks=np.array([i*0.1 for i in range(1,10)])
    subticks=list(subticks)+list(subticks*10)+list(subticks*100)
    '''
    在yscale的symlog中，它并不始终都是按照对数去画的，靠近原点的部分（-linthresh,linthresh）实际上是线性坐标。而linscale指定了该段线性坐标区间(0,linthresh)在图上的长度相对于对数坐标区间的长度。两者必须要设置得匹配，才能使得图像的坐标区间看起来是均匀的。python的default数值设置的是不匹配的。所以为了使坐标区间看起来均匀，得手动设置两个的值。
    '''
    plt.yscale('symlog',linthresh=0.1,linscale=1.0,subs=subticks)
    '''
    给一个建议，一般文献里画图时，fontsize都至少需要取到30以上，我一般会取到40（做PPT时也一样）。不然，坐标可能在文章中看不清楚。
    '''
    plt.xlabel(r"Voltage $(V)$", fontsize=35)
    plt.ylabel(r"Current $(\times 10^{-2}nA)$", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    plt.xlim([-20, 20])
    plt.ylim([-60, 60])
    plt.tick_params(axis='x', direction='in',width=2)
    plt.tick_params(axis='y', direction='in',width=2)
    plt.legend([r"$"+str(i)+"K$" for i in T], fontsize=15)
    '''
    将minor ticks的指向也设为'in'
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
    这儿
    是不是笔误了？传递的参数应该是x而不是T?
    HT_Plot(T, y, colormap='hsv',figname="HT_AVLog1.pdf")
    '''
    HT_Plot(x, y, colormap='hsv',figname="HT_AVLog1.pdf",T=T)
