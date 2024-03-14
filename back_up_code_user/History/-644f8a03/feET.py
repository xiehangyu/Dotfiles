import os
import sys
names=[]
for assoc in ['2','4','8']:
    for name in ['RandomRP','NMRURP','LRURP','LIPRP']:
        filename=name+"_"+assoc+'8'
        names.append(filename)

command="zip -9 -r 8kb.zip "+' '.join(names)
os.system(command)