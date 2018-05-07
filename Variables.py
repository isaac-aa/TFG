#--------------------------------------------
#		Variables.py
# This module allocate and store the different
# arrays of conservative and primite variables,
# as well as other auxiliary arrays that depend
# on the mesh shape
#
#--------------------------------------------
import numpy as np

import Parameters as par
import Grid
from encodings.bz2_codec import bz2_decode

# Conservative variables
rho = np.ones(Grid.z.shape)
momentumZ = np.zeros(Grid.z.shape)
momentumY = np.zeros(Grid.z.shape)
momentumX = np.ndarray([])
if par.twoHalfD==True:
   momentumX = np.zeros(Grid.z.shape)
energy = np.ones(Grid.z.shape)

# Conservative variables at last time step
lastrho = np.ones(Grid.z.shape)
lastmomentumZ = np.zeros(Grid.z.shape)
lastmomentumY = np.zeros(Grid.z.shape)
lastmomentumX = np.ndarray([])
if par.twoHalfD==True:
   lastmomentumX = np.zeros(Grid.z.shape)
lastenergy = np.ones(Grid.z.shape)

# Source arrays
momentumSourceZ = np.zeros(Grid.z.shape)
momentumSourceY = np.zeros(Grid.z.shape)
energySource = np.zeros(Grid.z.shape)

# Primitive variables
T = np.ones(Grid.z.shape)
P = np.ones(Grid.z.shape)
vZ = np.zeros(Grid.z.shape)
vY = np.zeros(Grid.z.shape)
vX = np.ndarray([])
if par.twoHalfD==True:
   vX = np.zeros(Grid.z.shape)

# Kappa of thermal diffusion
kappa = np.zeros(Grid.z.shape)

# Arrays for tridiagonal matrix
diag = np.zeros(Grid.z.shape)
lower = np.zeros(Grid.z.shape)
upper = np.zeros(Grid.z.shape)
rhs = np.zeros(Grid.z.shape)

RadiativeLoss = np.zeros(Grid.z.shape)
logT_table, Lamda_table = np.loadtxt('dere_etal_table.dat', usecols=(0,1), unpack=True)
logLamda_table = np.log10(Lamda_table)


Bz = np.ndarray([])
By = np.ndarray([])
Bx = np.ndarray([])
lastBz = np.ndarray([])
lastBy = np.ndarray([])
lastBx = np.ndarray([])
if par.MHD == True:
   Bz = np.zeros(Grid.z.shape)
   By = np.zeros(Grid.z.shape)
   Bx = np.zeros(Grid.z.shape)
   lastBz = np.zeros(Grid.z.shape)
   lastBy = np.zeros(Grid.z.shape)
   lastBx = np.zeros(Grid.z.shape)
   

def updateLastVars():
   lastrho[:] = rho[:]
   lastmomentumZ[:] = momentumZ[:]
   if par.dim==2:
      lastmomentumY[:] = momentumY[:]
   lastenergy[:] = energy[:]
   if par.MHD:
      lastBz[:] = Bz[:]
      lastBy[:] = By[:]
      lastBx[:] = Bx[:]
      lastmomentumX[:] = momentumX[:]

