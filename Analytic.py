#--------------------------------------------
#		Analytic.py
# This module defines some analytic functions
# to compare them to numerical solutions. They
# take as arguments the same as the initial
# condition plus the time
#
#--------------------------------------------

import numpy as np

import Grid
import Parameters as par
import Settings as sets
import Variables as var

def SoundWaves(args, t):
  rho0 = args[0]
  A = args[1]
  p0 = args[2]
  Nwav = args[3]

  #K = Nwav*2*np.pi/(Grid.zf-Grid.z0)
  c_s = np.sqrt(par.gamma*p0/rho0)
  phas = Nwav * 2. * np.pi *(Grid.z-par.z0)/(par.zf-par.z0)
  
  # v = w/k -> w = v*k
  w = c_s*Nwav*2.*np.pi/(par.zf-par.z0)
  
  rhoAna = rho0*(1. + A*np.cos(phas - w*t) )
  vAna =  c_s*A*np.cos(phas - w*t)
  PAna = p0*(1. + par.gamma*A*np.cos(phas - w*t) )
  TAna = PAna/(rhoAna*(par.gamma-1.)) 
 
  return rhoAna, vAna, PAna, TAna
  

def Isothermal(args, t):
  p0 = args[0]
  rho0 = args[1]
  PAna = p0*np.exp(-Grid.z/1.)
  rhoAna = rho0*np.exp(-Grid.z/1.)
  vAna = Grid.z*0.
  
  return rhoAna, vAna, PAna  
  
def GaussianThermal(args, t):
  T0 = args[0]
  z0 = args[1]
  width = args[2]
  rho0 = args[3]
  
  t_new = t + width/par.ct
  exp = np.sqrt(width/(par.ct*t_new))*np.exp(-(Grid.z-z0)*(Grid.z-z0)/( 4*par.ct*t_new) ) 
  
  rhoAna = rho0*np.ones(Grid.z.shape)
  e = (1.+T0*exp)*par.cv
  
  vAna = Grid.z*0.
  PAna = rhoAna*(par.gamma-1.)*e
  TAna = (1.+T0*exp) 

  return rhoAna, vAna, PAna, TAna
   
   
   
   
   
   
   
   
   
   
   
  
