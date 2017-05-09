#--------------------------------------------
#		Settings.py
# This module defines the different functions
# and numerical schemes used for the simulation
#
#--------------------------------------------

import numpy as np
import Advance
import BoundaryConditions
import InitialConditions
import ChangeOfVar
import Parameters as par

print 'Loading Settings..'


InitialCondition = InitialConditions.IsothermalEq
argsIC = [1., 1.]
#InitialCondition = InitialConditions.SoundWaves
#argsIC = [1.0, 0.001, 1.0, 4.]

ChangeOfVar.ConvertToPrim()


#BoundaryConditionL = BoundaryConditions.Periodic
#argsL = []
#BoundaryConditionR = BoundaryConditions.Periodic
#argsR = []

BoundaryConditionL = BoundaryConditions.FixedRhoP
argsL = ['L',1.,1.]
BoundaryConditionR = BoundaryConditions.Wall
argsR = ['R']


Scheme = Advance.FirstGen
#Advance.AllocateFirstGen() #Ponerlo en el main, despues de crear la malla


