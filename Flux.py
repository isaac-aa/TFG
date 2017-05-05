import numpy as np


print 'Loading Flux..'


# ------------------ FLUX ----------------------

def ComputeFlux(rho, momentum, energy, v, P):
   massFlux = momentum
   momentumFlux = rho*v*v + P
   energyFlux = (energy + P)*v
   return massFlux, momentumFlux, energyFlux











