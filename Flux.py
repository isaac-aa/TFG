import numpy as np
import ChangeOfVar
from Variables import * 

print 'Loading Flux..'


# ------------------ FLUX ----------------------

def ComputeFlux(rho, momentum, energy, v, P):
   massFlux = momentum
   momentumFlux = rho*v*v + P
   energyFlux = (energy + P)*v
   return massFlux, momentumFlux, energyFlux



def LaxFriedichs(dt, dz, rho, momentum, momentumSource, energy, energySource):
   lamda = dt/dz
   print 'LaxFriedichs'
   # Compute change of variables
   v = momentum/rho
   
   e = energy/rho - 0.5*v*v              #rho*e = E - 0.5*rho*v*v
   
   #T = e/Cv          
   P = rho*(gamma-1.)*e

   massFlux, momentumFlux, energyFlux = ComputeFlux(rho, momentum, energy, v, P)

   rho[1:-1] = 0.5*(rho[2:]+rho[:-2]) - 0.5*lamda*(massFlux[2:] - massFlux[:-2]) 
   momentum[1:-1] = 0.5*(momentum[2:]+momentum[:-2]) - 0.5*lamda*(momentumFlux[2:] - momentumFlux[:-2]) + dt*momentumSource[1:-1]
   energy[1:-1] = 0.5*(energy[2:]+energy[:-2]) - 0.5*lamda*(energyFlux[2:] - energyFlux[:-2]) + dt*energySource[1:-1]

   return rho, momentum, energy



def LaxFriedichsHalf(dt, dz, rho, momentum, momentumSource, energy, energySource):
   lamda = dt/dz

   # Compute change of variables
   v = momentum/rho
   
   e = energy/rho - 0.5*v*v              #rho*e = E - 0.5*rho*v*v
   
   #T = e/Cv          
   P = rho*(gamma-1.)*e

   massFlux, momentumFlux, energyFlux = ComputeFlux(rho, momentum, energy, v, P)

   rho[:-1] = 0.5*(rho[1:]+rho[:-1]) - 0.5*lamda*(massFlux[1:] - massFlux[:-1]) 
   momentum[:-1] = 0.5*(momentum[1:]+momentum[:-1]) - 0.5*lamda*(momentumFlux[1:] - momentumFlux[:-1]) + 0.5*dt*momentumSource[:-1]
   energy[:-1] = 0.5*(energy[1:]+energy[:-1]) - 0.5*lamda*(energyFlux[1:] - energyFlux[:-1]) + 0.5*dt*energySource[:-1]


   return rho, momentum, energy



def FirstGen(dt, dz, rho, momentum, momentumSource, energy, energySource):
   lamda = dt/dz

   rhoHalf = rho
   momentumHalf = momentum
   energyHalf = energy

   rhoHalf, momentumHalf, energyHalf = LaxFriedichsHalf(dt, dz, rhoHalf, momentumHalf, momentumSource, energyHalf, energySource) 

   # Compute change of variables
   vHalf = momentum/rho
   
   eHalf = energy/rho - 0.5*v*v              #rho*e = E - 0.5*rho*v*v
   
   #T = e/Cv          
   PHalf = rho*(gamma-1.)*eHalf


   massFluxHalf, momentumFluxHalf, energyFluxHalf = ComputeFlux(rhoHalf, momentumHalf, energyHalf, vHalf, PHalf)


   #rho[1:-1] = rho[1:-1] - lamda*( massFluxHalf[1:-1] - massFluxHalf[:-2] )
   #momentum[1:-1] = momentum[1:-1] - lamda*( momentumFluxHalf[1:-1] - momentumFluxHalf[:-2] ) + dt*momentumSource[1:-1]
   #energy[1:-1] = energy[1:-1] - lamda*( energyFluxHalf[1:-1] - energyFluxHalf[:-2] ) + dt*energySource[1:-1]
   rho[1:-1] = rho[1:-1] + lamda*( massFluxHalf[1:-1] - massFluxHalf[:-2] )
   momentum[1:-1] = momentum[1:-1] + lamda*( momentumFluxHalf[1:-1] - momentumFluxHalf[:-2] ) + dt*momentumSource[1:-1]
   energy[1:-1] = energy[1:-1] + lamda*( energyFluxHalf[1:-1] - energyFluxHalf[:-2] ) + dt*energySource[1:-1]

   
   return rho, momentum, energy









