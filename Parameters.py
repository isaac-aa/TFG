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
tf = 250.
it = 0
max_it = np.inf
cfl_set = 0.8
cfl = cfl_set

# ------------------ PLOT SETUP ---------------

# Default configuration
rhoAxis = []
vAxis = []
PAxis = []
TAxis = []
logScale = [False, False, False, False]

# Soundwaves
SoundSpeedLine = False
SoundSpeedAnalytic = False
#rhoAxis = [0.995,1.005]
#vAxis = [-.005, .005]
#PAxis = [0.995, 1.005]

# Isothermal eq
IsothermalAnalytic = False
#rhoAxis = [0., 1.1]
#vAxis = [-1., 0.1]
#PAxis = [0., 1.1]

# Transition region
rhoAxis = [1e-15, 5.7181381e-13]
vAxis = [-5e3, 5e3]
PAxis = [0., .5]
TAxis = [9000., 1100000.]
logScale = [True, False, False, True]

save_rate = 100


# ------------------ SOURCE TERMS -------------

# Gravity
IsComputingSource = True
IsThereGravity = True
g = -1.

# Momentum Damping
MomentumDamping = True
DampingPercent = 0.05

# ------------------ EQUATION OF STATE --------


gamma = 5./3.
R = 1.
Na = 6.022e23
molarMass = 1.6605e-24 #grams
mu = 1.
cv = 1.

def Computecv():
  global cv
  cv = R/(mu*(gamma-1.)*molarMass*Na)









