import sys
import csv
filename=sys.argv[1]
f=open(filename,'r')
lines=f.readlines()
f.close()
datas=[]
dataline=['','daxpy','unrolled daxpy','daxsbxpxy','unrolled daxsbxpxy','stencil','unrolled stencil']
datas.append(dataline)
key1=0
for line in lines:
    dataline1=['执行时常（CPU周期数）']
    if "simTicks" in line:
        print(key1)
        if key1>1:
            line=line.split()
            print(line)
            dataline1.append(line[1])
datas.append(dataline1)
with open(filename.replace(".txt",".csv"), 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(datas)

    