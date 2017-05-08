#--------------------------------------------
#		Flux.py
# This module computes the conservative flux
# given a state (rho, momentum, energy + v, P)
#
#--------------------------------------------

import numpy as np


print 'Loading Flux..'



def ComputeFlux(rho, momentum, energy, v, P):
   massFlux = momentum
   momentumFlux = rho*v*v + P
   energyFlux = (energy + P)*v
   return massFlux, momentumFlux, energyFlux











