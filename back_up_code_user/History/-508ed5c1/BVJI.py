import os
datas=[]
datas2=[]
datas3=[]
dataline=['','NMRURP','RandomRP','LRURP','LIPRP']
datas.append(dataline)
for assoc in ['2','4','8']:
    dataline0=[assoc]
    for name in ['nmru','randomru','lru','lip']:
        filename=name+"_"+assoc+".txt"
        f=open(filename+"/stats.txt",'r')
        lines=f.readlines()
        f.close()
        for line in lines:
            if "system.cpu.dcache.overall_misses::total" in line:
                line=line.split()
                dataline0.append(line[1])