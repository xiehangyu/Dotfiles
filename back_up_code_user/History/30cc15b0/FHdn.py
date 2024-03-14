import sys
import csv
filename=sys.argv[1]
f=open(filename,'r')
lines=f.readlines()
f.close()
datas=[]
dataline=['','daxpy','unrolled daxpy','daxsbxpxy','unrolled daxsbxpxy','stencil','unrolled stencil']
datas.append(dataline)
key0=0
dataline0=['CPI']
key1=0
dataline1=['执行时常（CPU周期数）']
key2=0
dataline2=['指令条数']
for line in lines:
    if "cpu.cpi" in line:
        key0=key0+1
        if key0>1 and key0<8:
            line=line.split()
            print(line)
            dataline0.append(line[1])
for line in lines:
    if "simTicks" in line:
        key1=key1+1
        if key1>1 and key1<8:
            line=line.split()
            print(line)
            dataline1.append(line[1])
for line in lines:
    if "cpu.numInsts" in line:
        key2=key2+1
        if key2>1 and key2<8:
            line=line.split()
            print(line)
            dataline2.append(line[1])
datas.append(dataline1)
with open(filename.replace(".txt",".csv"), 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(datas)

    