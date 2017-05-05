import numpy as np

def ComputeDT(dt, dt_max, dz, rho, momentum, cfl_set):
   cfl = dt*np.max(abs(momentum/rho))/dz
   if np.max(abs(momentum))!=0:
     dt = cfl_set*dz/np.max(abs(momentum/rho))
 
   if dt>=dt_max:
      dt = dt_max

   return cfl, dt


