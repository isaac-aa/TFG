import numpy as np
import Advance
import BoundaryConditions
import InitialConditions
import ChangeOfVar
import Variables as var

print 'Loading Settings..'

InitialConditions.IsothermalEq(1.)
InitialConditions.SoundWaves(1.0, 0.001, 1.0, 4.)


ChangeOfVar.ConvertToPrim()


BoundaryConditionL = BoundaryConditions.Periodic
argsL = []
BoundaryConditionR = BoundaryConditions.Periodic
argsR = []


"""
BoundaryConditionL = BoundaryConditions.FixedRhoP
argsL = ['L',1.,1.]
BoundaryConditionR = BoundaryConditions.Wall
argsR = ['R']
"""


Scheme = Advance.FirstGen



