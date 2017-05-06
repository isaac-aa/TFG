import numpy as np
import Variables as var

print "Loading TimeStep.."

def ComputeDT():

   
   c_s = np.sqrt(var.gamma*var.P/var.rho)
   vchar1 = np.max(abs(var.v + c_s))
   vchar2 = np.max(abs(var.v - c_s))
   
   var.cfl = var.dt*np.max( [vchar1, vchar2] )/var.dz

   if np.max(abs(var.v))!=0:
     var.dt = var.cfl_set*var.dz/np.max( [vchar1, vchar2] )
 
   if var.dt>=var.dt_max:
      var.dt = var.dt_max

