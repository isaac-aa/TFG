import numpy as np
from Variables import * 

print 'Loading Flux..'


def LaxFriedichs(dt, dz, rho, massFlux, momentum, momentumFlux, momentumSource, energy, energyFlux, energySource):
   lamda = dt/dz[1:-1]

   T = (energy - 0.5*momentum*momentum/rho)/(rho*Cv) 

   rho[1:-1] = 0.5*(rho[2:]+rho[:-2]) - 0.5*lamda*(massFlux[2:] - massFlux[:-2]) 
   
   momentum[1:-1] = 0.5*(momentum[2:]+momentum[:-2]) - 0.5*lamda*(momentumFlux[2:] - momentumFlux[:-2]) + dt*momentumSource[1:-1]

   energy[1:-1] = 0.5*(energy[2:]+energy[:-2]) - 0.5*lamda*(energyFlux[2:] - energyFlux[:-2]) + dt*energySource[1:-1]

   return rho, momentum, energy
