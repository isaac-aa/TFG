#--------------------------------------------
#      read.py
# Module for reading the data given by the simulation,
# Works for both 1D and 2D cases and can be easily configurable
#
#--------------------------------------------



import numpy as np
import matplotlib.pyplot as plt
import glob

case = 'OrszagTangVortex'

mainfolder = '../RESULTS/'+case+'/'
folder_dat = '../RESULTS/'+case+'/RESULTS_DAT/'
folder_im =  '../RESULTS/'+case+'/RESULTS/'

files = sorted(glob.glob(folder_dat+'*.npy'))
files = [ifile.split('/')[-1] for ifile in files]

print files

dim = 2

time = np.loadtxt(mainfolder+'log.txt', usecols=(2,))

#import Grid
#import Parameters as par
#Grid.Uniform1DGrid(par.Nz, par.z0, par.zf)
#import Analytic

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
      data = np.load(folder_dat+files[i])
      print i, ' out of ', len(files)

      gamma = 5./3.
      R = 1.
      mu = 1.
      

      vZ = 2.*data[3][1:-1,:-1]/(data[2][1:-1,:-1] + data[2][1:-1,1:])
      vY = 2.*data[4][:-1,1:-1]/(data[2][:-1,1:-1] + data[2][1:,1:-1])
      vX = data[5]/data[2]
      momZ = 0.5*(data[3][1:-1,:-2]+data[3][1:-1,1:-1])
      momY = 0.5*(data[4][:-2,1:-1]+data[4][1:-1,1:-1])
      P = np.zeros(data[0].shape)
      P[1:-1,1:-1] = (data[5][1:-1,1:-1] - 0.5*(momZ*momZ+momY*momY)/data[2][1:-1,1:-1])*(gamma-1.)
      P[1:-1,   0] = (data[5][1:-1,0] - 0.5*(momZ[:,0]*momZ[:,0]+data[4][1:-1, 0]*data[4][1:-1, 0])/data[2][1:-1,0])*(gamma-1.) # L
      P[1:-1,  -1] = (data[5][1:-1,-1] - 0.5*(momZ[:,-1]*momZ[:,-1]+data[4][1:-1,-2]*data[4][1:-1,-2])/data[2][1:-1,-1])*(gamma-1.) # R
      P[0,   1:-1] = (data[5][0,1:-1] - 0.5*(data[3][0,1:-1]*data[3][0,1:-1]+momY[0,:]*momY[0,:])/data[2][0,1:-1])*(gamma-1.) # B
      P[-1,  1:-1] = (data[5][-1,1:-1] - 0.5*(data[3][-2,1:-1]*data[3][-2,1:-1]+momY[-1,:]*momY[-1,:])/data[2][-1,1:-1])*(gamma-1.) # T
   
      P += (- 0.5*(data[-1]*data[-1] + data[-2]*data[-2] + data[-3]*data[-3]) - 0.5*data[-4]*data[-4]/data[2] )* (gamma-1.)
   
      T = P*mu/(data[2]*R)
            
      plt.pcolormesh(data[0], data[1], data[2], cmap='jet')
      
#     #plt.pcolormesh(data[0], data[1], data[2]-1., vmin=-3e-3, vmax=3e-3)
#     plt.plot(data[1][:,100], data[2][:,100]-1., 'r')
#     plt.plot(data[1][:,100], P[:,100]-1., 'g')
#     plt.plot(data[1][:,100], data[3][:,100]/np.sqrt(gamma), 'b')
#     plt.plot(data[1][:,100], data[-3][:,100]-1., 'm')
#     #plt.plot(data[0][3,1:-1], vY[3,:]/np.sqrt(gamma), 'b--')
#     #plt.plot(data[0][3,:], vX[3,:]/np.sqrt(gamma), 'b:')
#     v_ms = np.sqrt(gamma + 1.)
#     plt.axvline(np.sqrt(gamma)*t) 
#     plt.axvline(v_ms*t, c='m') 

      v_A = 1.
#     plt.ylim(-2e-3,2e-3)
#     plt.plot(data[1][:,100], data[-3][:,100], 'k')
#     plt.plot(data[1][:,100], data[-2][:,100], 'k--')
#     plt.plot(data[1][:,100], data[-1][:,100], 'k:')
#     plt.plot(0.5*(data[1][:-1,100]+data[1][1:,100]), vY[:,100]/v_A, 'r:')

      #Q = plt.quiver(data[0], data[1], data[-3], data[-2])
      plt.title(time[i])
      #plt.colorbar()


      #plt.plot(T[5,:])
      #plt.semilogy()
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
      plt.savefig(folder_im+'%d.png'%i)
      plt.clf()
      #plt.close(fig)
      
      
      
      
      
      
      
      
      
      
      
      
