import numpy as np
from Variables import * 

print 'Loading SourceTerm..'

def computeGravSource(rho, momentum):
   momentumGravSource = rho*g
   energyGravSource = momentum*g

   return momentumGravSource, energyGravSource

def ComputeSource(rho, momentum):
  energyS = np.zeros(z.shape)
  momentumS = np.zeros(z.shape) 
  if IsThereGravity:
     momentumG, energyG = computeGravSource(rho, momentum)
     momentumS += momentumG
     energyS += energyG
  return momentumS, energyS
