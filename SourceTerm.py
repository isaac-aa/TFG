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
import Variables as var

print 'Loading SourceTerm..'

def computeGravSource():
   momentumGravSource = var.rho*par.g
   energyGravSource = var.momentum*par.g
 
   return momentumGravSource, energyGravSource

def ComputeSource():
  if par.IsThereGravity:
     momentumG, energyG = computeGravSource()
     var.momentum += par.dt*momentumG
     var.energySource += par.dt*energyG

