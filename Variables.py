import numpy as np


# ------------------ MESH ---------------------
N = 256
z0 = 0.
zf = 1.

dz = (zf-z0)/N   #np.ones(z.shape)*(zf-z0)/N
z = np.linspace(z0-dz/2., zf+dz/2., N+2)

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

cfl_set = 0.001

SoundSpeedLine = True

dt_max = 0.01
dt = 1e-9
tt = 0.
tf = .75
it = 0
max_it = 2 #50000
save_rate = 5




