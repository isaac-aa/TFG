import numpy as np


# ------------------ MESH ---------------------
N = 20000
z0 = 0.
zf = 1.

z = np.linspace(z0, zf, N)
dz = np.ones(z.shape)*(zf-z0)/N

rho = np.ones(z.shape)
momentum = np.ones(z.shape)
energy = np.ones(z.shape)

massFlux = rho
momentumFlux = momentum
energyFlux = energy

momentumSource = np.zeros(z.shape)
energySource = np.zeros(z.shape)

T = np.ones(z.shape)



# ------------------ GENERAL CONFIGURATION ----

IsComputingSource = True
IsThereGravity = True
g = -1.

Cv = 1.
R = 1.

cfl_set = 0.05



dt_max = 0.00001
dt = dt_max
tt = 0.
tf = 1.
it = 0
max_it = 50000




