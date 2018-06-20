#--------------------------------------------
#      Settings.py
# This module defines the different functions
# and numerical schemes used for the simulation
#
#--------------------------------------------

import Advance
import BoundaryConditions as BC
import InitialConditions
import ChangeOfVar
import LagrangeTracer
import Parameters as par
import Derivative
import Grid

print 'Loading Settings..'

InitialCondition = InitialConditions.OrszagTangVortex
argsIC = []



ChangeOfVar.ConvertToPrim()


BoundaryConditionT = BC.Periodic(BC.BoundaryCondition('T'), BC.BoundaryCondition('B'))
BoundaryConditionB = None
BoundaryConditionT.setup()

BoundaryConditionR = BC.Periodic(BC.BoundaryCondition('R'), BC.BoundaryCondition('L'))
BoundaryConditionL = None
BoundaryConditionR.setup()


Tracers = None 

Derivate = Derivative.CentralDer2D()

Scheme = Advance.RK3()
Scheme.setup(Derivate)
