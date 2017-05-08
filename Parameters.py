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
tf = 0.8
it = 0
max_it = 800
cfl_set = 0.8
cfl = cfl_set

# ------------------ PLOT SETUP ---------------

save_rate = 1
SoundSpeedLine = False
SoundSpeedAnalytic = True

# ------------------ SOURCE TERMS -------------

IsComputingSource = False
IsThereGravity = False
g = -1.

# ------------------ EQUATION OF STATE --------

Cv = 1.
gamma = 5./3.
R = 1.










