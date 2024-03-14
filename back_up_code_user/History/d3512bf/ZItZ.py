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
        ����Ϊ�˰�ȫ��������Ǽ���fp.close()�ɡ�
        python��C��һ��������ļ����رյĻ�����ôԭ��д���ļ������ݾͲ�����µ�Ӳ���ϡ������Ļ���������������նˣ�
        �ͻᵼ��ԭ��д�������Ҳ��ʧ�ˡ���֮ǰ�ͷ������������⡣
        '''
        fp.close()
    return x, y

def HT_Plot(x,y,colormap='hsv',figname="HT_AVLog1.pdf",T=None):
    plt.figure(figsize=(10,8))
    '''
    ����һ������ȡcolormap
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
    ��log��ͼʱ�����黹�Ǽ���minor ticks�����subticks���ɡ�������˼�ǣ�������ʾ����1E-1,1E0,1E1�������⣬������ָ�����������ϣ�����2��3��4��������ʾ���ߡ��Ҹо���㻹��ͦ���õģ����԰�����log��ͼ�ϸ����׵ض�����
    '''
    plt.minorticks_on()
    '''
    ����minor ticks(���subticks)��λ�á�
    '''
    subticks=np.array([i*0.1 for i in range(1,10)])
    subticks=list(subticks)+list(subticks*10)+list(subticks*100)
    '''
    ��yscale�е�symlog��ʵ����һֱ�ǰ��ն������껭ͼ�ġ��ڿ���ԭ��ĸ���(-linthresh,linthresh)�ǰ����������껭ͼ�ġ�linscale������ָ���������������(0,linthresh)��ͼ�ϵĳ�������ڶ������굥�����䣨��������10-100���ı��������linthresh��linscale���õ�ֵ��ƥ�䣬 �ͻ�ʹ���������俴�������Ǿ��ȷֲ��ġ���python��default value�У������ߵ������ǲ�ƥ��ġ�Ϊ��ʹ���������俴�������ȷֲ�����Ҫ�ֶ��������ߵ�ֵ
    '''
    plt.yscale('symlog',linthresh=0.1,linscale=1.0,subs=subticks)
    '''
   ��һ��Ŀ���������ͼ�У��ҽ���fontsize����ȡ��30�ɣ���һ��ȡ��40����PPTʱҲһ��������Ȼ���ܻ��������￴�������
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
    �ҽ��黹��legend�ɣ��������Ժܷ���ÿ���ÿһ�����߶�Ӧʲô�¶ȡ�
    '''
    plt.legend([r"$"+str(i)+"K$" for i in T], fontsize=20)
    '''
    ʹminor ticks��ָ��ҲΪ'in'
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
    ����ǲ��Ǳ����ˣ��Ҹо����ﴫ��ȥ�Ĳ���Ӧ����x,������T?
    HT_Plot(T, y, colormap='hsv',figname="HT_AVLog1.pdf")
    '''
    HT_Plot(x, y, colormap='hsv',figname="HT_AVLog1.pdf",T=T)
    '''
    ����һ��T����������дlegend��
    '''
