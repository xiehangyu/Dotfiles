import numpy as np
import matplotlib.pyplot as plt
x1=[6,7,8,9,10]
x2=[6,7,8,9,10,11,12]
ybase=[909,2242,18722,193070,4693856]
yAXA=[110,520,3633,29879,260436,3041078,24948709]
YAXAblock=[111,521,3829,34050,293240,2267471,19925876]
YAXAproved=[43,213,1792,14130,163387,18475026]
plt.rcParams.update({'font.size': 22})
plt.figure(figsize=(12,8))
ybase=np.array(ybase)
yAXA=np.array(yAXA)
YAXAblock=np.array(YAXAblock)
YAXAproved=np.array(YAXAproved)
ybase=np.log2(ybase)
yAXA=np.log2(yAXA)
YAXAblock=np.log2(YAXAblock)
YAXAproved=np.log2(YAXAproved)
plt.plot(x1,ybase,'r')
plt.plot(x2,yAXA,'b')
plt.plot(x2,YAXAblock,'g')
plt.plot(x2,YAXAproved,'y')
plt.legend(['baseline','basic AVX',"AVX Block",'AVX Improved','','','',''])
plt.scatter(x1,ybase,c='r')
plt.scatter(x2,yAXA,c='b')
plt.scatter(x2,YAXAblock,c='g')
plt.scatter(x2,YAXAproved,c='y')
plt.xlabel(r'$\log_2$N')
plt.ylabel(r'$\log_2$Time')