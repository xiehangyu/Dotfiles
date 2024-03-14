import os
datas=[]
datas2=[]
datas3=[]
dataline=['','NMRURP','RandomRP','LRURP','LIPRP']
datas.append(dataline)
for assoc in ['2','4','8']:
    dataline0=[assoc]
    for name in ['nmru','randomru','lru','lip']:
        filename=