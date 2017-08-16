#--------------------------------------------
#      Parameters.py
# In this module all the different user defined
# parameters are set.
#
#--------------------------------------------

import numpy as np

print "Loading Variables.."

# ------------------ MESH ---------------------
Nz = 5000
z0 = 1.351040753526313007e+08
zf = 9.184177197366108894e+08


Ny = 1
y0 = -5e7
yf = 5e7

dim = 1

staggered = True

# ----------------- TIME SETUP ----------------

dt_max = 0.1
dt = dt_max
tt = 0.
tf = np.inf #50000.
it = 0
max_it = np.inf
cfl_set = .9
cfl = cfl_set

# ------------------ PLOT SETUP ---------------


FolderName = 'RESULTS/AnalyticalLossEq_HydrostaticP_NoDamping_SmoothLambda'
save_rate = 10000


# ------------------ SOURCE TERMS -------------


IsComputingSource = True

# Gravity
IsThereGravity = True
GravityMode = 'Constant' 
FreeFallLine = False
g = 1.

# Momentum Damping
MomentumDamping = False
DampingMode = 'MaximumVelocity'
DampingMultiplier = 0.01
DampingMaxVel = 1000000

# Thermal Diffusion
ThermalDiffusion = True
SpitzerDiffusion = True
ImplicitConduction = True
ct = 9e-12 * 1e5 # Value at E.Priest "Solar Magnetohydrodynamics" * mks to cgs factor
f_cfl = .5  #0.5 Largest stable number
DiffusionPercent = 1. #1e-2



# Radiative losses
RadiativeLoss = True
RadiationPercent = 1.

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








