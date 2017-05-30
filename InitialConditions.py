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
  var.momentum = var.momentum*0.
  var.energy = P0/(par.gamma-1.)


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
  var.momentum = vIn*var.rho

  pIn = p0*(1. + par.gamma*A*np.cos(phas) )
  var.energy = pIn/(par.gamma-1.) + 0.5*var.rho*vIn*vIn


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
  var.momentum = var.momentum*0.
  var.rho = np.exp(logrho)
  var.energy = var.rho*par.cv*var.T


def LogTGravityProfile(args):
  logTA = np.log(args[0])
  logTB = np.log(args[1])
  rhoA = args[2]

  par.g = -27360.00
  par.R = 8.3144598e+07
  par.mu = 1.113202

  logT = logTA + (logTB-logTA)/(Grid.z[-1]-Grid.z[0]) * (Grid.z-Grid.z[0])
  var.T = np.exp(logT)
 
  logp = np.zeros(Grid.z.shape)
  logp[0] = np.log(rhoA*par.R*args[0]/par.mu)

  B = par.mu*par.g/par.R
  A = par.ct
  Lambda = 0.

  i=1
  # Compute pressure to make it stationary under gravity
  while i<len(Grid.z):

    logp[i] = logp[i-1] + Grid.dz*B/var.T[i]   

    i+=1  

  print np.exp(logp)
  par.Computecv()
  var.rho = np.exp(logp)*par.mu/(par.R * var.T)
  var.momentum = var.momentum*0.
  var.energy = var.rho*par.cv*var.T


def GaussianTemperature(args):
  T0 = args[0]
  z0 = args[1]
  width = args[2]
  rho0 = args[3]
  
  
  exp = np.exp(-(Grid.z-z0)*(Grid.z-z0)/(4*width))
  var.rho = rho0*np.ones(Grid.z.shape)
  var.momentum = var.momentum*0.
  var.energy = (1.+T0*exp)*var.rho*par.cv
  if par.SpitzerDiffusion:
    var.kappa = par.ct*np.ones(Grid.z.shape)*(1.+T0*exp)**(5./2.)
  else: 
    var.kappa = par.ct*np.ones(Grid.z.shape) #*exp #np.ones(Grid.z.shape)
  
  
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
   var.momentum = var.momentum*0.
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
   
   p_data, rho_data = np.loadtxt(FileName, skiprows=2, usecols=(1,2), unpack=True)

   var.rho = rho_data*rho_ref
   var.energy = p_data /(par.gamma-1.)
   var.momentum = var.momentum*0.
   var.T = var.energy/(var.rho*par.cv)
   var.kappa = par.ct*(var.T)**(5./2.)
   #print var.kappa
   
def RestartFromFile(args):
   print "Restarting simulation..."
   files = glob.glob(args[0]+'/*.dat')
   files.sort()
   
   last_it = files[-1]

   f = open(last_it)
   refs = f.readline().split()
   mu = float(refs[1])
   g = abs(float(refs[2]))
   R = float(refs[3])
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
