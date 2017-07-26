import numpy as np
import matplotlib.pyplot as plt
import glob

folder = '../RESULTS/Star/RESULTS_DAT/'

files = glob.glob(folder+'*.npy')

for i in files:
  data = np.load(i)
  T = (5./3.-1.)*(data[5] - 0.5*(data[3]*data[3] + data[4]*data[4])/data[2])
  plt.imshow(data[2])  #data[5]- 0.5*(data[3]*data[3] + data[4]*data[4])/data[2])
  #plt.plot(data[0][10,:], data[2][10,:]) 
  #plt.ylim(1+2e-3, 1-2e-3) 
 
  plt.savefig(i+'.png')

  plt.cla()
