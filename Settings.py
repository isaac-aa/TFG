import numpy as np
import Flux
import BoundaryConditions
import InitialConditions
from Variables import * 

print 'Loading Settings..'

rho, momentum, energy = InitialConditions.IsothermalEq(z, rho, momentum, energy)

BoundaryConditionL = BoundaryConditions.FixedRho
argsL = ['L',1.]
BoundaryConditionR = BoundaryConditions.Wall
argsR = ['R']


FluxScheme = Flux.LaxFriedichs



