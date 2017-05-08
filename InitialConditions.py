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

# ------------------ INITIAL CONDITIONS --------

print 'Loading InitialConditions..'


def IsothermalEq(args):
  p0 = args[0]
  P0 = p0*np.exp(-Grid.z/1.)
  var.rho = 1.*np.exp(-Grid.z/1.)
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
  
  
  
