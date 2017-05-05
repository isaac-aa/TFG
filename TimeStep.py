import numpy as np
import Variables as var

print "Loading TimeStep.."

def ComputeDT():
   var.cfl = var.dt*np.max(abs(var.v))/var.dz
   if np.max(abs(var.v))!=0:
     var.dt = var.cfl_set*var.dz/np.max(abs(var.v))
 
   if var.dt>=var.dt_max:
      var.dt = var.dt_max

