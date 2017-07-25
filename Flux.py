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




def ComputeFlux2D(rho, momentumZ, momentumY, energy, vZ, vY, P):
   massFluxZ = momentumZ
   massFluxY = momentumY
   momentumZFluxZ = rho*vZ*vZ + P
   momentumZFluxY = rho*vZ*vY
   momentumYFluxZ = rho*vZ*vY
   momentumYFluxY = rho*vY*vY + P
   energyFluxZ = (energy + P)*vZ
   energyFluxY = (energy + P)*vY
   return massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY







