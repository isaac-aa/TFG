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
#InitialCondition = InitialConditions.ReadICFromFilePressure
#argsIC = ['Extras/ThermalEq_IC_2.dat']
#InitialCondition = InitialConditions.RestartFromFile
#argsIC = ['RESULTS/AnalyticalLossEq_HydrostaticP/RESULTS_DAT']
#InitialCondition = InitialConditions.ReadICFromFile
#argsIC = ['hydrostatic_equilibrium_2.dat']
#InitialCondition = InitialConditions.GaussianTemperature
#argsIC = [2., 0.5, 0.0005, 1.]
#InitialCondition = InitialConditions.IsothermalEq
#argsIC = [1., 1.]
#InitialCondition = InitialConditions.SoundWaves2D
#argsIC = [1.0, 0.001, 1.0, 4.]
#InitialCondition = InitialConditions.PureMagnetosonic2D
#argsIC = [1.0, 0.001, 1.0, 4., 1.]
#InitialCondition = InitialConditions.PureAlfven2D
#argsIC = [1.0, 0.001, 1.0, 1.0,4.]
#InitialCondition = InitialConditions.KelvinHelmholtz2D
#argsIC = [1e-2, 1e-2, .3, -.3, 0.5, 0.25, 0.05, 50]
#InitialCondition = InitialConditions.IsothermalAtm2D
#argsIC = [1., 1.]
#InitialCondition = InitialConditions.RayleighTaylorIns
#argsIC = [1.,2., 0.5, 10., 0.2, .5]



ChangeOfVar.ConvertToPrim()


"""
BoundaryConditionR = BC.BCComposite('R')
BoundaryConditionL = BC.BCComposite('L')

BoundaryConditionR.setup(BC.ZeroDer, (BC.AntiSym,BC.AntiSym), BC.FixedT, (BC.ZeroDer,BC.ZeroDer,BC.ZeroDer), [None, None,None, 1e6])
BoundaryConditionL.setup(BC.ZeroDer, (BC.AntiSym,BC.AntiSym), BC.FixedT, (BC.ZeroDer,BC.ZeroDer,BC.ZeroDer), [None, None, None, 1e5])


BoundaryConditionT = BC.BCComposite('T')
BoundaryConditionB = BC.BCComposite('B')

BoundaryConditionT.setup(BC.ZeroDer, (BC.AntiSym,BC.AntiSym), BC.ZeroDer, (BC.ZeroDer,BC.ZeroDer,BC.ZeroDer),[None, None, None, None])
BoundaryConditionB.setup(BC.ZeroDer, (BC.AntiSym,BC.AntiSym), BC.ZeroDer, (BC.ZeroDer,BC.ZeroDer,BC.ZeroDer),[None, None, None, None])
"""


BoundaryConditionT = BC.Periodic(BC.BoundaryCondition('T'), BC.BoundaryCondition('B'))
BoundaryConditionB = None
BoundaryConditionT.setup()

BoundaryConditionR = BC.Periodic(BC.BoundaryCondition('R'), BC.BoundaryCondition('L'))
BoundaryConditionL = None
BoundaryConditionR.setup()


Tracers = None #(LagrangeTracer.Tracer(par.z0 + 500*Grid.dz, 0), LagrangeTracer.Tracer(par.z0 + 1000*Grid.dz, 0), LagrangeTracer.Tracer(par.z0 + 1500*Grid.dz, 0), LagrangeTracer.Tracer(par.z0 + 2000*Grid.dz, 0) , LagrangeTracer.Tracer(par.z0 + 2500*Grid.dz, 0) , LagrangeTracer.Tracer(par.z0 + 3000*Grid.dz, 0) , LagrangeTracer.Tracer(par.z0 + 3500*Grid.dz, 0))

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
#BoundaryConditionL = BoundaryConditions.ConservativeMassHydrostaticP #WallExpRhoHydrostaticP #WallFixedRhoFixedT #FixedT #WallSecondRhoFixedT
#argsL = ['L', 1e5, 3.7677286546362324e-14]
#BoundaryConditionR = BoundaryConditions.FixedT
#argsR = ['R', 1e6, 2.52373286e-15]

Derivate = Derivative.CentralDer2D()

Scheme = Advance.RK3()
Scheme.setup(Derivate)
