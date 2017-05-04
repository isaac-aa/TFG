import numpy as np

# ------------------ INITIAL CONDITIONS --------

print 'Loading InitialConditions..'


def IsothermalEq(z, rho, momentum, energy):
  return 1.*np.exp(-z/1.), momentum*0., np.ones(z.shape)
