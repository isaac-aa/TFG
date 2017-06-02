#--------------------------------------------
#		Variables.py
# This module allocate and store the different
# arrays of conservative and primite variables,
# as well as other auxiliary arrays that depend
# on the mesh shape
#
#--------------------------------------------
import numpy as np

import Grid

# Conservative variables
rho = np.ones(Grid.z.shape)
momentum = np.zeros(Grid.z.shape)
energy = np.ones(Grid.z.shape)

# Conservative variables at last time step
lastrho = np.ones(Grid.z.shape)
lastmomentum = np.ones(Grid.z.shape)
lastenergy = np.ones(Grid.z.shape)

# Source arrays
momentumSource = np.zeros(Grid.z.shape)
energySource = np.zeros(Grid.z.shape)

# Primitive variables
T = np.ones(Grid.z.shape)
P = np.ones(Grid.z.shape)
v = np.ones(Grid.z.shape)

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




