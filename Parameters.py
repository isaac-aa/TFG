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
tf = 1.
it = 0
max_it = np.inf
cfl_set = 0.9
cfl = cfl_set

# ------------------ PLOT SETUP ---------------

FolderName = 'RESULTS_TESTCASES/SoundWaves_FirstGen' #'RESULTS/AnalyticalLossEq_HydrostaticP_ConsRho_SymV_NoDamping'
save_rate = 1

SaveToFile = True
SaveToFileRatio = 1

PlotFile = False
FileToPlot = 'Extras/ThermalLossesEq_ana.dat'

PlotCharacteristics = False

# Default configuration
rhoAxis = []
vAxis = []
PAxis = []
TAxis = []

ylabels = [r'$\frac{g}{cm^3}$', r'$\frac{cm}{s}$', r'$\frac{dyn}{cm^2}$', r'$K$']
xlabel = r'$cm$'
logScale = [False, False, False, False]

# Soundwaves
SoundSpeedLine = False
SoundSpeedAnalytic = True
rhoAxis = [0.995,1.005]
vAxis = [-.005, .005]
PAxis = [0.995, 1.005]
TAxis = [1.495, 1.505]

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
#rhoAxis = [6e-16, 5e-12]
#vAxis = [-10e5, 10e5]
#PAxis = [0., .5]
#TAxis = [2e3, 1.1e6]
#logScale = [True, False, False, True]


# ------------------ SOURCE TERMS -------------


IsComputingSource = False

# Gravity
IsThereGravity = False
FreeFallLine = False
g = -1.

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
Na = 1 #6.022e23
molarMass = 1 #1.6605e-24 #grams
mu = 1.
cv = 1.

def Computecv():
  global cv
  cv = R/(mu*(gamma-1.)*molarMass*Na)









