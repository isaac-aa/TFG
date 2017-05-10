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
   momentumGravSource = 0.5*(var.rho+var.lastrho)*par.g
   energyGravSource = 0.5*(var.momentum+var.lastmomentum)*par.g
   return momentumGravSource, energyGravSource


def computeMomentumDamping():
   DampingVel = par.DampingPercent*var.v
   
   momentumDampingSource = -DampingVel*var.rho
   energyDampingSource = -0.5*DampingVel*DampingVel*var.rho
   return momentumDampingSource, energyDampingSource

def ComputeSource():
  if par.IsThereGravity:
     momentumG, energyG = computeGravSource()
     var.momentum += par.dt*momentumG
     var.energySource += par.dt*energyG
  if par.MomentumDamping:
     momentumDamping, energyDamping = computeMomentumDamping()
     var.momentum += par.dt*momentumDamping
     var.energySource += par.dt*energyDamping
