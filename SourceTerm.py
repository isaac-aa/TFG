#--------------------------------------------
#		Source.py
# This module computes the momentum and energy
# source term that should be added to the 
# conservative part of the equations. This can be
# further improve using Operator Splitting (TODO)
#
#--------------------------------------------

import numpy as np
import Parameters as par

print 'Loading SourceTerm..'

def computeGravSource():
   momentumGravSource = var.rho*par.g
   energyGravSource = var.momentum*par.g
 
   return momentumGravSource, energyGravSource

def ComputeSource():
  if par.IsThereGravity:
     momentumG, energyG = computeGravSource()
     var.momentumSource = momentumG
     var.energySource = energyG

