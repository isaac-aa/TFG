#--------------------------------------------
#		TimeStep.py
# This module in on charge of computing the 
# time step each iteration
#
#--------------------------------------------

import numpy as np
import Parameters as par
import Grid
import Variables as var

print "Loading TimeStep.."

def ComputeDT():

   
   c_s = np.sqrt(par.gamma*var.P/var.rho)
   vchar1 = np.max(abs(var.v + c_s))
   vchar2 = np.max(abs(var.v - c_s))
   
   par.cfl = par.dt*np.max( [vchar1, vchar2] )/Grid.dz

   par.dt = par.cfl_set*Grid.dz/np.max( [vchar1, vchar2] )
 
