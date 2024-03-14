import sys
filename=sys.argv[1]
f=open(filename,'r')
lines=f.readlines()
f.close()
f=open(filename.replace(".txt",".csv"),'w')
datas=[]
dataline=['','daxpy','unrolled daxpy','daxsbxpxy','unrolled daxsbxpxy','stencil','unrolled stencil']
datas.append(dataline)
for line in lines:
    key1=0
    dataline1=['执行时常（CPU周期数）']
    if "simTicks" in line:
        print(line)
        key1+=1
        if key>1:
            line=line.split()
            f.write(line[1]+",")
f.write("\n")
f.close()
    