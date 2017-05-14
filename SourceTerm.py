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
import Grid

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

def computeTemperatureDiffusion():
   #ct = 9e-12 * 1e5 # Value at E.Priest "Solar Magnetohydrodynamics" * mks to cgs factor
   if par.SpitzerDiffusion:
      var.kappa = par.ct*var.T**5./2.
   
   #tau = Grid.dz*Grid.dz*var.rho[-2]*par.cv/var.kappa
   #print 'tau: %3.e'%np.min(tau)
   
   Dkappa = 0.5*(var.kappa[2:]-var.kappa[:-2])
   EnergyDiff = (var.kappa[1:-1]+Dkappa)*var.T[2:] - 2.*var.kappa[1:-1]*var.T[1:-1] + (var.kappa[1:-1]-Dkappa)*var.T[:-2]

   return EnergyDiff/(Grid.dz*Grid.dz)



def computeRadiativeLosses():
   logT = np.log10(var.T)
   Q = np.zeros(Grid.z.shape)  #This can be allocated at settings/variables
   for i in range(len(logT)):   #Not sure if this will be faster than masked arrays...
     if logT[i] <= 4.97:
       Q[i] = 1.09e-31*var.T[i]*var.T[i]
     elif logT[i] < 5.67:
       Q[i] = 8.87e-17/var.T[i]
     elif logT[i] < 6.18:
       Q[i] = 1.90e-22
     elif logT[i] < 6.55:
       Q[i] = 3.53e-13*var.T[i]**(-3./2.)
     elif logT[i] < 6.90:
       Q[i] = 3.56e-25*var.T[i]**(1./3.)
     else:
       Q[i] = 5.49e-16/var.T[i]
   
   return var.rho * Q 


def ComputeSource():
  if par.IsThereGravity:
     momentumG, energyG = computeGravSource()
     var.momentum += par.dt*momentumG
     var.energy += par.dt*energyG
  if par.MomentumDamping:
     momentumDamping, energyDamping = computeMomentumDamping()
     var.momentum += momentumDamping
     var.energy += energyDamping
  if par.ThermalDiffusion:
     ThermalDiff = computeTemperatureDiffusion()
     var.energy[1:-1] += par.dt*ThermalDiff
  if par.RadiativeLoss:
     RadLoss = computeRadiativeLosses()
     var.energy -= par.dt*RadLoss
     
     
     
     
     
     
     
     
