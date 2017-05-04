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

massFlux = np.zeros(z.shape)
momentumFlux = np.zeros(z.shape)
energyFlux = np.zeros(z.shape)

momentumSource = np.zeros(z.shape)
energySource = np.zeros(z.shape)

T = np.ones(z.shape)
P = np.ones(z.shape)
v = np.ones(z.shape)

# ------------------ GENERAL CONFIGURATION ----

IsComputingSource = False
IsThereGravity = False
g = -1.

Cv = 1.
gamma = 5./3.
R = 1.

cfl_set = 0.05



dt_max = 0.00001
dt = 1e-9
tt = 0.
tf = 1.
it = 0
max_it = 50000




