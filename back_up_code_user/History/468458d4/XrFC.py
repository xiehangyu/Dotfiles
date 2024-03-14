import os
for i in range(3,12):
    for j in range(3,8):
        blocksize=(1<<j)
        nsize=i
        print("Now for nsize={}, blocksize={}".format(nsize,blocksize))
        command="./Lab5CPU.out "+str(nsize)+" "+str(blocksize)
        os.system(command)