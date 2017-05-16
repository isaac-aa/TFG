import numpy as np

import SourceTerm
import Parameters as par
import Variables as var
import Grid


def Thermal():
  DT_z = (var.T[1:]-var.T[:-1])/Grid.dz
  L = np.abs(var.T[:-1]/DT_z)

  Chi = var.kappa/(var.rho*par.cv)
  DiffusionTime = (5./2. + 1.)*L*L/Chi[:-1]
  
  return DiffusionTime, L


def Radiative():
  Lr = SourceTerm.computeRadiativeLosses()
  
  RadiationTime = par.cv*var.T*var.rho/Lr

  return RadiationTime


def Dynamic():
  c_s = np.sqrt(par.gamma*var.P/var.rho)
  SoundTime = Grid.dz/c_s
  
  return SoundTime

