import numpy as np

print "Loading Variables.."

# ------------------ MESH ---------------------
N = 100
z0 = 0.
zf = 1.
dz = (zf-z0)/N

dz_p = (zf-z0)/(N+2)   #np.ones(z.shape)*(zf-z0)/N
z = np.linspace(z0, zf, N+1)
z = np.append([z0-dz], z)   #Todo esto no es que sea muy optimo...
z = np.append(z, [zf+dz])
print z, dz

rho = np.ones(z.shape)
momentum = np.ones(z.shape)
energy = np.ones(z.shape)

"""
massFlux = np.zeros(z.shape)
momentumFlux = np.zeros(z.shape)
energyFlux = np.zeros(z.shape)
"""

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

cfl_set = 0.97
cfl = cfl_set

SoundSpeedLine = True

dt_max = 0.1
dt = dt_max
tt = 0.
tf = 0.8
it = 0
max_it = 800
save_rate = 1




