import numpy as np
from Variables import *
# ------------------ INITIAL CONDITIONS --------

print 'Loading InitialConditions..'


def IsothermalEq(z, rho, momentum, energy, p0):
  P0 = p0*np.exp(-z/1.)
  rho0 = 1.*np.exp(-z/1.)
  e = P0/((gamma-1.)*rho)
  return rho0, momentum*0., e


def SoundWaves(z, rho, momentun, energy, rho0, rho1, p0, N):
  K = N*2*np.pi/(z[1]-z[-2])
  rhoIn = rho0 + rho1*np.cos(K*z)
  
  c_s = np.sqrt(gamma*p0/rho0)
  vIn =  c_s*rho1*np.cos(K*z)

  pIn = p0 + gamma*rho1*np.cos(K*z)
  EIn = pIn/(gamma-1.) + 0.5*rhoIn*vIn*vIn

  return rhoIn, vIn*rhoIn, EIn
