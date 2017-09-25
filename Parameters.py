#--------------------------------------------
#      Parameters.py
# In this module all the different user defined
# parameters are set.
#
#--------------------------------------------

import numpy as np

print "Loading Variables.."

# ------------------ MESH ---------------------
Nz = 288
z0 = 0
zf = 2.*np.pi


Ny = 288
y0 = 0
yf = 2.*np.pi

dim = 2
twoHalfD = True
MHD = True

staggered = True

# ----------------- TIME SETUP ----------------

dt_max = 0.1
dt = dt_max
tt = 0.
tf = 3. #np.inf
it = 0
max_it = np.inf
cfl_set = .75
cfl = cfl_set

# ------------------ PLOT SETUP ---------------


FolderName = 'RESULTS/OrszagTangVortex'
save_rate = 1


# ------------------ SOURCE TERMS -------------


IsComputingSource = False

# Gravity
IsThereGravity = False
GravityMode = 'Constant' 
FreeFallLine = False
g = 1.

# Momentum Damping
MomentumDamping = False
DampingMode = 'MaximumVelocity'
DampingMultiplier = 0.01
DampingMaxVel = 1000000

# Thermal Diffusion
ThermalDiffusion = False
SpitzerDiffusion = False
ImplicitConduction = False
ct = 9e-12 * 1e5 # Value at E.Priest "Solar Magnetohydrodynamics" * mks to cgs factor
f_cfl = .5  #0.5 Largest stable number
DiffusionPercent = 1. #1e-2



# Radiative losses
RadiativeLoss = False
RadiationPercent = 1.

# ------------------ EQUATION OF STATE --------


gamma = 5./3.
R = 1.
Na = 1. #6.022e23
molarMass = 1. #1.6605e-24 #grams
mu = 1.
cv = 1.

def Computecv():
   global cv
   cv = R/(mu*(gamma-1.)*molarMass*Na)








