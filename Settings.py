import numpy as np
import Flux
import BoundaryConditions
import InitialConditions
import ChangeOfVar
from Variables import * 

print 'Loading Settings..'

#rho, momentum, energy = InitialConditions.IsothermalEq(z, rho, momentum, energy, 1.)
rho, momentum, energy = InitialConditions.SoundWaves(z, rho, momentum, energy, 1.0, 0.01, 1.0, 3.)


v,T,P = ChangeOfVar.ConvertToPrim(rho, momentum, energy)


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


FluxScheme = Flux.FirstGen



