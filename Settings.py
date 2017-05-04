import numpy as np
import Flux
import BoundaryConditions
import InitialConditions
from Variables import * 

print 'Loading Settings..'

#rho, momentum, energy = InitialConditions.IsothermalEq(z, rho, momentum, energy, 1.)
rho, momentum, energy = InitialConditions.SoundWaves(z, rho, momentum, energy, 1.0, 0.01, 1.0, 3.)


# Compute change of variables
v = momentum/rho
 
e = energy/rho - 0.5*v*v              #rho*e = E - 0.5*rho*v*v
   
T = e/Cv          
P = rho*(gamma-1.)*e


BoundaryConditionL = BoundaryConditions.Periodic
argsL = []
BoundaryConditionR = BoundaryConditions.Periodic
argsR = []


"""
BoundaryConditionL = BoundaryConditions.FixedRho
argsL = ['L',1.]
BoundaryConditionR = BoundaryConditions.Wall
argsR = ['R']
"""


FluxScheme = Flux.LaxFriedichs



