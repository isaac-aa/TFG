#--------------------------------------------
#		Parameters.py
# In this module all the different user defined
# parameters are set.
#
#--------------------------------------------

import numpy as np

print "Loading Variables.."

# ------------------ MESH ---------------------
N = 5000
z0 = 1.349473185696847439e+08 #0.
zf = 9.185744765195573568e+08 #1.

# ----------------- TIME SETUP ----------------

dt_max = 0.1
dt = dt_max
tt = 0.
tf = np.inf #500.
it = 0
max_it = np.inf
cfl_set = 0.1
cfl = cfl_set

# ------------------ PLOT SETUP ---------------

save_rate = 1000

SaveToFile = True
SaveToFileRatio = 10

PlotFile = True
FileToPlot = 'Extras/ThermalEq.dat'

PlotCharacteristics = True

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

# Thermal Diffusion
ThermalAnalytic = False
#rhoAxis = [0., 2.]
#vAxis = [-1., 1.]
#PAxis = [0., 2.]
#TAxis = [0., 3.]

# Transition region
SoundSpeedProfile = False
#rhoAxis = [1e-15, 5e-12]
#vAxis = [-1e5, 1e5]
#PAxis = [0., 2.]
#TAxis = [1e3, 1.8e6]
#logScale = [True, False, False, True]

# Thermal Equilibrium
rhoAxis = [1e-15, 5e-12]
vAxis = [-1e5, 1e5]
PAxis = [0., .5]
TAxis = [8e4, 1.5e6]
logScale = [True, False, False, True]


# ------------------ SOURCE TERMS -------------


IsComputingSource = True

# Gravity
IsThereGravity = True
g = -1.

# Momentum Damping
MomentumDamping = True
DampingMultiplier = 1.

# Thermal Diffusion
ThermalDiffusion = True
SpitzerDiffusion = True
ImplicitConduction = True
ct = 9e-12 * 1e5 # Value at E.Priest "Solar Magnetohydrodynamics" * mks to cgs factor
f_cfl = .5  #0.5 Largest stable number
DiffusionPercent = 1. #1e-2



# Radiative losses
RadiativeLoss = False
RadiationPercent = 1e-2

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









