import sys
filename=sys.argv[1]
f=open(filename,'r')
lines=f.readlines()
f.close()
f=open(filename.replace("./txt","./csv"),'w')