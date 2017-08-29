#--------------------------------------------
#      read.py
# Module for reading the data given by the simulation,
# Works for both 1D and 2D cases and can be easily configurable
#
#--------------------------------------------



import numpy as np
import matplotlib.pyplot as plt
import glob

folder = '../RESULTS/2D_Staggered_STR_Implicit/'

files = glob.glob(folder+'RESULTS_DAT/*.npy')

dim = 2

import Grid
import Parameters as par
Grid.Uniform1DGrid(par.Nz, par.z0, par.zf)
import Analytic

if dim == 1:
   for i in range(len(files)):
      data = np.load(files[i])
      t = float(files[i][-30:-8])  #0.00000000000000000000
      print i, ' out of ', len(files), t
      
      fig, axs = plt.subplots(ncols=1, nrows=4)      
      fig.set_size_inches(18.5, 10.5)

      axs[0].plot(data[0], data[2])
      axs[0].set_title(r'$\rho$')
      axs[0].semilogy()
      axs[0].set_ylim(1e-15, 1e-12) 
      
      
      axs[1].plot(data[0], data[3])
      axs[1].set_title(r'$v_z$')
      mod = np.max(np.abs(data[3]))
      axs[1].set_ylim(-mod, mod)
      
      R = 8.3144598e+07
      mu = 1.113202
      axs[2].plot(data[0], data[2]*data[5]*R/mu)
      axs[2].set_title(r'$P$')
      axs[2].set_ylim(0, .5)
      
      axs[3].plot(data[0], data[5])
      axs[3].set_title(r'$T$')
      axs[3].semilogy()
      axs[3].set_ylim(5e3, 1.1e6)

      """
      TracerZ, TracervZ = np.loadtxt(files[i][:-8] + '.tcr', unpack = True)
      for j in range(len(TracerZ)):
         for ax in axs:
            ax.axvline(x=TracerZ[j])
      """
      plt.tight_layout()
      plt.savefig(files[i]+'.png')

      plt.clf()
      plt.close()

if dim == 2:
   for i in range(len(files)):
      data = np.load(files[i])
      print i, ' out of ', len(files)
      
      #plt.imshow(data[2][:-1,1:-1], aspect='auto')
      #plt.title(r'$T$')
      #plt.colorbar()

      plt.plot(data[5][5,:])
      plt.semilogy()
      """
      fig, axs = plt.subplots(ncols=3, nrows=2)

      fig.set_size_inches(18.5, 10.5)


      img = axs[0,0].pcolor(data[0], data[1], data[2])
      axs[0,0].set_title(r'$\rho$')
      plt.colorbar(img, ax=axs[0,0])

      img = axs[0,1].pcolor(data[0], data[1], data[5]*data[2])
      axs[0,1].set_title(r'P')
      plt.colorbar(img, ax=axs[0,1])


      img = axs[0,2].pcolor(data[0], data[1], data[5])
      plt.colorbar(img, ax=axs[0,2])
      axs[0,2].set_title(r'T')
  
      img = axs[1,0].pcolor(data[0], data[1], data[3])
      axs[1,0].set_title(r'$v_z$')
      plt.colorbar(img, ax=axs[1,0])


      img = axs[1,1].pcolor(data[0], data[1], data[4])
      axs[1,1].set_title(r'$v_y$')
      plt.colorbar(img, ax=axs[1,1])

  
      axs[1,2].quiver(data[0], data[1], data[3], data[4])
      
      #plt.imshow(data[2])  #data[5]- 0.5*(data[3]*data[3] + data[4]*data[4])/data[2])
      #plt.plot(data[0][10,:], data[2][10,:]) 
      #plt.ylim(1+2e-3, 1-2e-3) 
      """
      plt.tight_layout()
      plt.savefig(files[i]+'.png')
      plt.clf()
      #plt.close(fig)
      
      
      
      
      
      
      
      
      
      
      
      