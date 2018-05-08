#--------------------------------------------
#		InitialConditions.py
# This module defines the Initial Conditions 
# for the simulation. They are selected at
# Settings.py, and a "args" argument is passed,
# which is a list of the neccesary user defined
# parameters.
#
#--------------------------------------------


import numpy as np
import shutil

import Grid
import Parameters as par
import Variables as var
import glob

# ------------------ INITIAL CONDITIONS --------

print 'Loading InitialConditions..'


def IsothermalEq(args):
   p0 = args[0]
   rho0 = args[1]
  
   P0 = p0*np.exp(-Grid.z/1.)
   var.rho = rho0*np.exp(-Grid.z/1.)
   var.momentumZ = var.momentumZ*0.
   var.energy = P0/(par.gamma-1.)


def ConstantFlow(args):
   T0 = args[0]
   rho0 = args[1]
   v = args[1]

   var.rho = rho0*np.ones(Grid.z.shape)
   var.momentumZ = var.rho*v
   var.energy = par.cv*T0*np.ones(Grid.z.shape) + 0.5*var.rho*v*v

def SoundWaves(args):
   rho0 = args[0]
   A = args[1]
   p0 = args[2]
   Nwav = args[3]

   #K = Nwav*2*np.pi/(Grid.zf-Grid.z0)
   phas = Nwav * 2. * np.pi *(Grid.z-par.z0)/(par.zf-par.z0)
  
   var.rho = rho0*(1. + A*np.cos(phas) )
  
   c_s = np.sqrt(par.gamma*p0/rho0)
   vIn =  c_s*A*np.cos(phas)
   var.momentumZ = vIn*var.rho

   pIn = p0*(1. + par.gamma*A*np.cos(phas) )
   var.energy = pIn/(par.gamma-1.) + 0.5*var.rho*vIn*vIn


def SoundWaves2D(args):
   rho0 = args[0]
   A = args[1]
   p0 = args[2]
   Nwav = args[3]

   phas = Nwav * 2. * np.pi *(Grid.z-par.z0)/(par.zf-par.z0)
  
   var.rho = rho0*(1. + A*np.cos(phas) )
  
   c_s = np.sqrt(par.gamma*p0/rho0)
   v_ms = np.sqrt(c_s**2 + var.Bx*var.Bx/var.rho)
   if par.staggered:
      phas2 = phas + np.pi*Grid.dz/(par.zf-par.z0)
      var.momentumZ = v_ms*A*np.cos(phas2)*np.ones_like(Grid.z)
   else:
      var.momentumZ = A*np.cos(phas)*np.ones_like(Grid.z)

   var.momentumY = var.momentumY*0.

   pIn = p0*(1. + par.gamma*A*np.cos(phas) )
   var.energy = pIn/(par.gamma-1.) + 0.5*var.momentumZ*var.momentumZ/var.rho


def PureMagnetosonic2D(args):
   rho0 = args[0]
   A = args[1]
   p0 = args[2]
   Nwav = args[3]
   B0 = args[4]

   phas = Nwav * 2. * np.pi *(Grid.y-par.y0)/(par.yf-par.y0)
  
   var.Bx = B0+B0*A*np.cos(phas)

   var.rho = rho0*(1. + A*np.cos(phas) )
  
   c_s = np.sqrt(par.gamma*p0/rho0)
   v_ms = np.sqrt(c_s**2 + var.Bx*var.Bx/var.rho)
   if par.staggered:
      phas2 = phas + np.pi*Grid.dz/(par.zf-par.z0)
      var.momentumY = v_ms*A*np.cos(phas2)*np.ones_like(Grid.z)
   else: 
      var.momentumY = v_ms*A*np.cos(phas)*np.ones_like(Grid.z)
   var.momentumZ = var.momentumZ*0.

   pIn = p0*(1. + par.gamma*A*np.cos(phas) )
   var.energy = pIn/(par.gamma-1.) + 0.5*var.momentumY*var.momentumY/var.rho + 0.5*var.Bx*var.Bx


def PureAlfven2D(args):
   rho0 = args[0]
   A = args[1]
   p0 = args[2]
   B0 = args[3]
   Nwav = args[4]

   phas = Nwav * 2. * np.pi *(Grid.y-par.y0)/(par.yf-par.y0)
   phas2 = phas + np.pi*Grid.dy/(par.yf-par.y0)

   var.rho = rho0*np.ones_like(Grid.z)
   var.By  = B0*np.ones_like(Grid.z)
   var.Bz  = B0*A*np.cos(phas)
   v_A = np.sqrt(B0*B0/rho0)
   var.momentumZ = -v_A*A*np.cos(phas)

   var.energy = p0/(par.gamma-1.) + 0.5*var.momentumZ*var.momentumZ/var.rho + 0.5*(var.Bz*var.Bz+var.Bx*var.Bx)


def KelvinHelmholtz2D(args):
   rho0 = args[0]
   p0 = args[1]
   vA = args[2]
   vB = args[3]
   center = args[4]
   width = args[5]
   amp = args[6]
   freq = args[7]
  
   mask_center = np.abs(Grid.y-center) < width+amp*np.sin(freq*Grid.z)/2.  #Mask for values INSIDE the A tube
   mask_outside = np.invert(mask_center)

   var.rho = np.ones(Grid.z.shape)*rho0

   var.momentumY = np.zeros(Grid.z.shape)
   var.momentumZ[mask_center] = vA*var.rho[mask_center]
   var.momentumZ[mask_outside] = vB*var.rho[mask_outside]

   var.energy = p0/(par.gamma-1.) + 0.5*var.momentumZ*var.momentumZ/var.rho 
   

def IsothermalAtm2D(args):
   p0 = args[0]
   rho0 = args[1]
  
   P0 = p0*np.exp(-Grid.z/1.)
   var.rho = rho0*np.exp(-Grid.z/1.)
   var.momentumY = var.momentumY*0.
   var.momentumZ = var.momentumZ*0.
   var.energy = P0/(par.gamma-1.)

   #pert = 1.*np.exp( -(Grid.z-0.5)*(Grid.z-0.5)/2. - (Grid.y-0.5)*(Grid.y-0.5)/2.  ) 
   #var.energy += pert



def RayleighTaylorIns(args):
   rho1 = args[0]
   rho2 = args[1]
   discZ = args[2]
   discP = args[3]
   amp = args[4]
   freq = args[5]


   # Mask for both domains
   thresholdZ = discZ - amp*np.sin(2.*np.pi*freq*Grid.y)
   upper = Grid.z>thresholdZ
   lower = np.invert(upper)

  
   # Static in all the domain
   var.momentumZ = var.momentumZ*0.
   var.momentumY = var.momentumY*0.

   # We set a constant density atmosphere at both domains
   var.rho[lower] = rho1
   var.rho[upper] = rho2

   P0 = -par.g*var.rho*(Grid.z-discZ) + discP
   var.energy = P0/(par.gamma-1.)


def BuildSystemMatrix(center, L, R,  T, B):
   N = Grid.z.shape[0]*Grid.z.shape[1]
   x_side = Grid.z.shape[0]    #Make sure that this is indeed correct for non-square meshes
   y_side = Grid.z.shape[1]


   nnz = np.zeros(5*N)  #Four neighbours for each cell plus the cell itself
   col = np.zeros(5*N, dtype=int)
   offset = np.zeros(N+1, dtype=int)
   i=1

   #center = 0. #-2*(1./(Grid.dz*Grid.dz) + 1./(Grid.dy*Grid.dy))
   #LR =   1./(dx)
   #TB = 0. #1./(Grid.dy*Grid.dy)

   for xi in range(Grid.z.shape[0]):
      for yj in range(Grid.z.shape[1]):
         index = (yj*y_side + xi)*5
         cell = yj*y_side + xi
         index = cell*5

         nnz[index] = center #Center cell
         col[index] = cell


         nnz[index + 1] = L #Left
         col[index + 1] = cell-1  #Left neighbor


         nnz[index + 2] = R #Right
         col[index + 2] = cell+1  #Right


         nnz[index + 3] = T #Top
         col[index + 3] = cell-x_side #Top


         nnz[index + 4] = B #Bottom
         col[index + 4] = cell+x_side


         offset[i] = offset[i-1]+5  #Improve performance-wise later on
         i+=1

   # Zero-derivative potential boundary condition
   for xi in range(Grid.z.shape[0]):
      index = (y_side*0 + xi)*5
      nnz[index + 3] = 0.
      col[index + 3] = N-(x_side-xi)
      nnz[index + 4] = -1.

      index = (y_side*(y_side-1) + xi)*5
      nnz[index + 4] = 0.
      col[index + 4] = xi
      nnz[index + 3] = -1.


   for yj in range(Grid.z.shape[1]):
      index = (y_side*yj + 0)*5
      cell = y_side*yj
      index = cell*5

      nnz[index + 1] = 0.
      col[index + 1] = cell + (x_side-1)
      nnz[index+2] = -1.

      index = (y_side*yj + (x_side-1))*5
      cell = y_side*yj + (x_side-1)
      nnz[index + 2] = 0.
      col[index + 2] = cell - (x_side-1)
      nnz[index+1] = -1



   return nnz, col, offset

def Star2D(args):
   """
   WIP
   """
   
   rho0 = args[0]
   rho1 = args[1]
   sigma = args[2]
   z0 = args[3]
   y0 = args[4]
   momY = args[5]
   p0 = 1.

   # Background density
   var.rho = np.ones(Grid.z.shape)*rho0

   # 'Planet'
   var.rho += rho1*np.exp(-((Grid.z-z0)*(Grid.z-z0) + (Grid.y-y0)*(Grid.y-y0) )/sigma  )
   #var.rho += rho1*np.exp(-((Grid.z-.2)*(Grid.z-.2) + (Grid.y-.2)*(Grid.y-.2) )/sigma  )

   # Static
   var.momentumZ = -momY*np.exp(-((Grid.z-z0)*(Grid.z-z0) + (Grid.y-y0)*(Grid.y-y0) )/sigma  )
   var.momemtumY = var.momentumY*0.

     
   fz = np.zeros(Grid.z.shape)
   fy = np.zeros(Grid.z.shape)

   import SourceTerm
   import scipy.sparse
   import scipy.sparse.linalg
   nnz, col, offset = SourceTerm.BuildSystemMatrix()
   N = Grid.z.shape[0]*Grid.z.shape[1]   #Repeated operation
   SystemMatrix = scipy.sparse.csr_matrix( (nnz, col, offset), shape=(N, N)  )
      
   #var.rho = 0.1*np.ones(Grid.z.shape)*np.exp( -( (Grid.z-0.4)*(Grid.z-0.4)+(Grid.y-0.4)*(Grid.y-0.4) )/0.005 )
   #var.rho += 0.1*np.exp( -( (Grid.z-0.7)*(Grid.z-0.7)+(Grid.y-0.7)*(Grid.y-0.7) )/0.005 )
      
   sol = scipy.sparse.linalg.spsolve(SystemMatrix, var.rho.flatten() )
   phi = sol.reshape(Grid.z.shape)
      
   fy[1:-1, 1:-1] = -var.rho[1:-1, 1:-1]*(phi[2:, 1:-1] - phi[:-2, 1:-1])/(2.*Grid.dz)
   fz[1:-1, 1:-1] = -var.rho[1:-1, 1:-1]*(phi[1:-1, 2:] - phi[1:-1, :-2])/(2.*Grid.dy)
   
   nnz, col, offset = BuildSystemMatrix(0., 1./(2.*Grid.dz), -1./(2.*Grid.dz), 0., 0.)
   N = Grid.z.shape[0]*Grid.z.shape[1]   #Repeated operation
   SystemMatrix = scipy.sparse.csr_matrix( (nnz, col, offset), shape=(N, N)  )   
   sol = scipy.sparse.linalg.spsolve(SystemMatrix, fz.flatten() )
   pZ = sol.reshape(Grid.z.shape)
 
   nnz, col, offset = BuildSystemMatrix(0., 0., 0., 1./(2.*Grid.dy), -1./(2.*Grid.dy) )
   N = Grid.z.shape[0]*Grid.z.shape[1]   #Repeated operation
   SystemMatrix = scipy.sparse.csr_matrix( (nnz, col, offset), shape=(N, N)  )   
   sol = scipy.sparse.linalg.spsolve(SystemMatrix, fy.flatten() )
   pY = sol.reshape(Grid.z.shape) 
   

   var.P = pZ*pY
   var.P += abs(np.min(var.P)) + 10.
   
   var.energy = var.P/(par.gamma-1.)
   
   import matplotlib.pyplot as plt
   
   plt.pcolor(Grid.z, Grid.y, var.P)
   plt.colorbar()
   plt.quiver(Grid.z, Grid.y, fz/var.rho, fy/var.rho)
   plt.show()





def LogTLogRhoProfile(args):
   logTA = np.log(args[0])
   logTB = np.log(args[1])
   rhoA = args[2]
   rhoB = args[3]

   par.g = -27360.00
   par.R = 8.3144598e+07
   par.mu = 1.113202

   logrho = np.log(rhoA) + (np.log(rhoB)-np.log(rhoA))/(Grid.z[-1]-Grid.z[0]) * (Grid.z-Grid.z[0])

   logT = logTA + (logTB-logTA)/(Grid.z[-1]-Grid.z[0]) * (Grid.z-Grid.z[0])
   var.T = np.exp(logT)

  
   """ 
   logp1 = np.zeros(var.P.shape)
   logp1[0] = logpA
  
   T_num = np.zeros(Grid.z.shape)
   T_num[0] = args[0]

   B = par.mu*par.g/par.R
   A = par.ct
   Lambda = 0.
 
   i=1
   while i<len(Grid.z):
      T7_2 = Grid.dz* Grid.dz * 7*Lambda/(2.*A) - T_num[i-2]**(7./2.) + 2*T_num[i-1]**(7./2.)
      T_num[i] = T7_2**(2./7.)

      logp[i] = logp[i-1] + Grid.dz*B/var.T[i]   
      logp1[i] = logp1[i-1] + Grid.dz*B/T_num[i]

      i+=1  
   """

   par.Computecv()
   var.momentumZ = var.momentumZ*0.
   var.rho = np.exp(logrho)
   var.energy = var.rho*par.cv*var.T


def LogTGravityProfile(args):
   logTA = np.log(args[0])
   logTB = np.log(args[1])
   rhoA = args[2]

   par.g = -27360.00
   par.R = 8.3144598e+07
   par.mu = 1.113202

   logT = logTA + (logTB-logTA)/(Grid.z[0,-1]-Grid.z[0,0]) * (Grid.z-Grid.z[0,0])
   var.T = np.exp(logT)
  
   logp = np.zeros(Grid.z.shape[1])
   logp[0] = np.log(rhoA*par.R*args[0]/par.mu)

   B = par.mu*par.g/par.R
   A = par.ct
   Lambda = 0.

   i=0
      
   # Compute pressure to make it stationary under gravity
   print logp.shape
   while i<par.Nz+2:

      logp[i] = logp[i-1] + Grid.dz*B/var.T[1,i]   

      i+=1  


   par.Computecv()
   
   var.rho = np.exp(logp)*par.mu/(par.R * var.T)
   var.momentumZ = var.momentumZ*0.
   var.energy = var.rho*par.cv*var.T
   if par.SpitzerDiffusion:
      var.kappa = par.ct*(var.T)**(5./2.)

def GaussianTemperature(args):
   T0 = args[0]
   z0 = args[1]
   width = args[2]
   rho0 = args[3]
  
  
   exp = np.exp(-(Grid.z-z0)*(Grid.z-z0)/(4*width))
   var.rho = rho0*np.ones(Grid.z.shape)
   var.momentumZ = var.momentumZ*0.
   var.energy = (1.+T0*exp)*var.rho*par.cv
   if par.SpitzerDiffusion:
      var.kappa = par.ct*np.ones(Grid.z.shape)*(1.+T0*exp)**(5./2.)
   else: 
      var.kappa = par.ct*np.ones(Grid.z.shape)*exp #np.ones(Grid.z.shape)
  
  
def ReadICFromFileTemperature(args):
   print "Loading IC from " + args[0]
   FileName = args[0]
   f = open(FileName)
   f.readline()     #Line containing number of elements
   refs = f.readline().split()
   T_ref = float(refs[1])
   rho_ref = float(refs[2])
   mu = float(refs[3])
   g = float(refs[4])
   R = float(refs[5])
   f.close()
   
   par.R = R
   par.mu = mu
   par.g = -g
   par.Computecv()
   #print par.R, par.mu, par.g, par.cv
   
   T_data, rho_data = np.loadtxt(FileName, skiprows=2, usecols=(1,2), unpack=True)

   var.rho = rho_data*rho_ref
   var.energy = var.rho*T_data*T_ref*R/mu /(par.gamma-1.)
   var.momentumZ = var.momentumZ*0.
   var.kappa = par.ct*(T_data*T_ref)**(5./2.)
   
   
def ReadICFromFilePressure(args):
   print "Loading IC from " + args[0]
   FileName = args[0]
   f = open(FileName)
   f.readline()     #Line containing number of elements
   refs = f.readline().split()
   p_ref = float(refs[1])
   rho_ref = float(refs[2])
   mu = float(refs[3])
   g = float(refs[4])
   R = float(refs[5])
   f.close()
   
   par.R = R
   par.mu = mu
   par.g = -g
   par.Computecv()
   #print par.R, par.mu, par.g, par.cv
   
   z_data, p_data, rho_data = np.loadtxt(FileName, skiprows=2, usecols=(0,1,2), unpack=True)

   var.rho =  np.interp(Grid.z, z_data, rho_data*rho_ref)
   var.energy = np.interp(Grid.z, z_data, p_data) /(par.gamma-1.)
   var.momentumZ = var.momentumZ*0.
   var.T = var.energy/(var.rho*par.cv)
   var.kappa = par.ct*(var.T)**(5./2.)
   #print var.kappa
   
   
   
def RestartFromFile(args):
   """ Deprecated
   files = glob.glob(args[0]+'/*.dat')
   files.sort()
  
   if len(args) == 2: #ie. there are two arguments, one of them is the filename 
      last_it = args[0] + '/' + args[1]
   else:
      last_it = files[-1]

   shutil.copy2(last_it, par.FolderName)

   print "Restarting simulation from " + last_it

   f = open(last_it)
   refs = f.readline().split()
   mu = float(refs[3])
   g = abs(float(refs[4]))
   R = float(refs[5])
   f.close()   
   print mu, R, g
   par.R = R
   par.mu = mu
   par.g = -g
   par.Computecv()
   
   rho_data, v_data, T_data = np.loadtxt(last_it, skiprows=1, usecols=(1,2,3), unpack=True)
  
   var.rho = rho_data
   var.energy = var.rho*T_data*par.cv + 0.5*var.rho*v_data*v_data
   var.momentum = v_data*var.rho
   var.kappa = par.ct*(T_data)**(5./2.)
   """
   
   
   files = glob.glob(args[0]+'/*.dat')
   files.sort()
  
   if len(args) == 2: #ie. there are two arguments, one of them is the filename 
      last_it = args[0] + '/' + args[1]
   else:
      last_it = files[-1]

   shutil.copy2(last_it, par.FolderName)   
   
   print "Restarting simulation from " + last_it 
      
   if par.dim==2:
      #We are assuming that the grid remains constant
      data = np.load(last_it)
      
      par.Computecv()
      
      var.rho = data[2]
      var.momentumZ = data[3]
      var.momentumY = data[4]
      var.energy = data[5]
      
      
def OrszagTangVortex(args):
   var.rho = np.ones(Grid.z.shape)*par.gamma*par.gamma
   var.P = np.ones(Grid.z.shape)*par.gamma
   
   var.momentumZ = -np.sin(Grid.y)*var.rho
   var.momentumY = np.sin(Grid.z)*var.rho
   var.momentumX = var.momentumX*0.
   
   var.Bz = -np.sin(Grid.y)
   var.By = np.sin(2*Grid.z)
   var.Bx = var.Bx*0.
   
   var.energy = var.P/(par.gamma-1.) + 0.5*var.rho*((-np.sin(Grid.y))*(-np.sin(Grid.y)) + np.sin(Grid.z)*np.sin(Grid.z)) + 0.5*(var.Bz*var.Bz + var.By*var.By)
   
   
   # Forced periodic conditions  (does not change anything...)
   """
   var.momentumZ[0,:] = var.momentumZ[-2,:]
   var.momentumZ[-1,:] = var.momentumZ[1,:]
   
   var.momentumY[:,0] = var.momentumY[:,-2]
   var.momentumY[:,-1] = var.momentumY[:,1]
   
   var.Bz[0,:] = var.Bz[-2,:]
   var.Bz[-1,:] = var.Bz[1,:]
   
   var.By[:,0] = var.By[:,-2]
   var.By[:,-1] = var.By[:,1]
   """
   #print var.momentumZ[1,:]/var.rho[1,:] - var.momentumZ[2,:]/var.rho[2,:]
   #print Grid.z[:,1]/(2.*np.pi), Grid.z[:,-2]/(2.*np.pi)
   #print var.momentumZ[1,:] - var.momentumZ[-2,:]
   
   
   
   
   
   
   
   
   
   
      
      
