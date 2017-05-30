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


#InitialCondition = InitialConditions.LogTGravityProfile
#argsIC = [1e5, 1e6, 3.7677286546362324e-14, 2.52373286e-15]
InitialCondition = InitialConditions.ReadICFromFilePressure
argsIC = ['Extras/ThermalEq_IC_2.dat']
#InitialCondition = InitialConditions.RestartFromFile
#argsIC = ['../Radiativ4/Parte1/RESULTS_DAT']
#InitialCondition = InitialConditions.ReadICFromFile
#argsIC = ['hydrostatic_equilibrium_2.dat']
#InitialCondition = InitialConditions.GaussianTemperature
#argsIC = [2., 0.5, 0.0005, 1.]
#InitialCondition = InitialConditions.IsothermalEq
#argsIC = [1., 1.]
#InitialCondition = InitialConditions.SoundWaves
#argsIC = [1.0, 0.001, 1.0, 4.]

ChangeOfVar.ConvertToPrim()


#BoundaryConditionL = BoundaryConditions.Periodic
#argsL = []
#BoundaryConditionR = BoundaryConditions.Periodic
#argsR = []

#BoundaryConditionL = BoundaryConditions.FixedRhoP
#argsL = ['L',1.,1.]
#BoundaryConditionR = BoundaryConditions.Wall
#argsR = ['R']

# Thermal diffusion
#BoundaryConditionR = BoundaryConditions.Wall
#argsR = ['R']
#BoundaryConditionL = BoundaryConditions.Wall
#argsL = ['L']


# Solar transition region
BoundaryConditionL = BoundaryConditions.WallFixedRhoFixedT #FixedT #WallSecondRhoFixedT
argsL = ['L', 1e5, 3.7677286546362324e-14]
BoundaryConditionR = BoundaryConditions.FixedT #WallFixedRhoFixedT
argsR = ['R', 1e6, 2.52373286e-15]


Scheme = Advance.FirstGen
#Advance.AllocateFirstGen() #Ponerlo en el main, despues de crear la malla


