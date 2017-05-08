import numpy as np

print "Loading Variables.."

# ------------------ MESH ---------------------
N = 100
z0 = 0.
zf = 1.
#dz = (zf-z0)/N
"""
#CreateMesh() #Desde el main
#Grid
dz_p = (zf-z0)/(N+2)   #np.ones(z.shape)*(zf-z0)/N
z = np.linspace(z0, zf, N+1)
z = np.append([z0-dz], z)   #Todo esto no es que sea muy optimo...
z = np.append(z, [zf+dz])
print z, dz
"""


# ------------------ GENERAL CONFIGURATION ----

IsComputingSource = False
IsThereGravity = False
g = -1.

Cv = 1.
gamma = 5./3.
R = 1.

cfl_set = 0.8
cfl = cfl_set

SoundSpeedLine = True

dt_max = 0.1
dt = dt_max
tt = 0.
tf = 0.8
it = 0
max_it = 800
save_rate = 1




