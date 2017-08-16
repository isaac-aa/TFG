#--------------------------------------------
#		Flux.py
# This module computes the conservative flux
# given a state (rho, momentum, energy + v, P)
#
#--------------------------------------------


print 'Loading Flux..'

import Parameters as par
import numpy as np
import Grid







def ComputeFlux(rho, momentum, energy):
   massFlux = momentum
   
   v2 = momentum*momentum/rho

   momentumFlux = (1.-0.5*(par.gamma-1.)) * v2 + energy*(par.gamma-1.) #rho*v*v + P
   
   energyFlux = (par.gamma*energy - (par.gamma-1.)*0.5*v2)*momentum/rho #(energy + P)*v

   return massFlux, momentumFlux, energyFlux




def ComputeFlux2D(rho, momentumZ, momentumY, energy):
   massFluxZ = momentumZ
   massFluxY = momentumY
   
   v2Z = momentumZ*momentumZ/rho
   v2Y = momentumY*momentumY/rho
   
   momentumZFluxZ = (1.-0.5*(par.gamma-1.)) * v2Z + energy*(par.gamma-1.) #rho*vZ*vZ + P
   momentumZFluxY = momentumZ*momentumY/rho #rho*vZ*vY
   
   momentumYFluxZ = momentumZFluxY          #rho*vZ*vY
   momentumYFluxY = (1.-0.5*(par.gamma-1.)) * v2Y + energy*(par.gamma-1.) #rho*vY*vY + P
   
   energyFluxZ = (par.gamma*energy - (par.gamma-1.)*0.5*(v2Z+v2Y)/rho)*momentumZ/rho  #(energy + P)*vZ
   energyFluxY = (par.gamma*energy - (par.gamma-1.)*0.5*(v2Z+v2Y)/rho)*momentumY/rho          #(energy + P)*vY
   
   return massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY







