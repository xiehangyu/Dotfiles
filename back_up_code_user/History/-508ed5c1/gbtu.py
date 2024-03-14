import os
import csv
datas=[]
datas2=[]
datas3=[]
dataline=['','NMRURP','RandomRP','LRURP','LIPRP']
datas.append(dataline)
datas2.append(dataline)
datas3.append(dataline)
for assoc in ['2','4','8']:
    dataline0=[assoc]
    dataline1=[assoc]
    dataline2=[assoc]
    for name in ['NMRURP','RandomRP','LRURP','LIPRP']:
        filename=name+"_"+assoc+'.txt'
        f=open(filename+"/stats.txt",'r')
        lines=f.readlines()
        f.close()
        for line in lines:
            if "dcache.demandMissRate" in line:
                line=line.split()
                dataline0.append(line[1])
                print(line[1])
                print("haha")
            if "cpu.cpi" in line:
                line=line.split()
                dataline1.append(line[1])
            if "simTicks" in line:
                line=line.split()
                dataline2.append(line[1])
    datas.append(dataline0)
    datas2.append(dataline1)
    datas3.append(dataline2)
with open("missrate8k.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(datas)
with open("cpi8k.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(datas2)
with open("simTicks8k.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(datas3)