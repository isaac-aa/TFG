import numpy as np
import Variables as var

print 'Loading SourceTerm..'

def computeGravSource():
   momentumGravSource = var.rho*var.g
   energyGravSource = var.momentum*var.g
 
   return momentumGravSource, energyGravSource

def ComputeSource():
  if var.IsThereGravity:
     momentumG, energyG = computeGravSource()
     var.momentumSource = momentumG
     var.energySource = energyG

