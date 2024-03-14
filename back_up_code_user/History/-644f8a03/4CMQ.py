import os
import sys
names=[]
for assoc in ['2','4','8']:
    for name in ['RandomRP','NMRURP','LRURP','LIPRP']:
        filename=name+"_"+assoc+'8'
        names.append(filename)
        newname=filename[0:-1]
        comand="cp -r "+filename+"  ./8kb/"+newname
        print(comand)
        os.system(comand)
