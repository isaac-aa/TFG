import numpy as np
import Variables as var

# ------------------ INITIAL CONDITIONS --------

print 'Loading InitialConditions..'


def IsothermalEq(p0):
  P0 = p0*np.exp(-var.z/1.)

  var.rho = 1.*np.exp(-var.z/1.)
  var.momentum = var.momentum*0.
  var.energy = P0/((var.gamma-1.)*var.rho)


def SoundWaves(rho0, A, p0, Nwav):
  #K = Nwav*2*np.pi/(var.zf-var.z0)
  phas = Nwav * 2. * np.pi *(var.z-var.z0)/(var.zf-var.z0)
  
  var.rho = rho0*(1. + A*np.cos(phas) )
  print var.rho[1]- var.rho[-2]
  print var.rho[0]- var.rho[-1]
  
  print var.rho[100]- var.rho[-101] 
  
  c_s = np.sqrt(var.gamma*p0/rho0)
  vIn =  c_s*A*np.cos(phas)
  var.momentum = vIn*var.rho

  pIn = p0*(1. + var.gamma*A*np.cos(phas) )
  var.energy = pIn/(var.gamma-1.) + 0.5*var.rho*vIn*vIn
