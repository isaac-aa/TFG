#--------------------------------------------
#		Parameters.py
# In this module all the different user defined
# parameters are set.
#
#--------------------------------------------

import numpy as np

print "Loading Variables.."

# ------------------ MESH ---------------------
N = 100
z0 = 0.
zf = 1.

# ----------------- TIME SETUP ----------------

dt_max = 0.1
dt = dt_max
tt = 0.
tf = 10.
it = 0
max_it = 150
cfl_set = 0.8
cfl = cfl_set

# ------------------ PLOT SETUP ---------------

# Soundwaves
SoundSpeedLine = False
SoundSpeedAnalytic = False
#rhoAxis = [0.995,1.005]
#vAxis = [-.005, .005]
#PAxis = [0.995, 1.005]

# Isothermal eq
IsothermalAnalytic = True
rhoAxis = [0., 1.1]
vAxis = [-1., 0.1]
PAxis = [0., 1.1]


save_rate = 1


# ------------------ SOURCE TERMS -------------

IsComputingSource = True
IsThereGravity = True
g = -1.

# ------------------ EQUATION OF STATE --------

Cv = 1.
gamma = 5./3.
R = 1.










