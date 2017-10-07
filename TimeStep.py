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
   vcharZ = np.max([np.max(abs(var.vZ + c_s)), np.max(abs(var.vZ - c_s))])
   
   if par.dim == 2:
      vcharY = np.max([np.max(abs(var.vY + c_s)), np.max(abs(var.vY - c_s))])
      par.dt = par.cfl_set/( vcharZ/Grid.dz + vcharY/Grid.dz  )
   else:
      par.dt = par.cfl_set*Grid.dz/vcharZ

   if par.ThermalDiffusion and not par.ImplicitConduction:
      dt_thermal = par.f_cfl*np.min(par.cv*var.rho*Grid.dz*Grid.dz/var.kappa) 
      #print 'Thermal: %.3e \t Sound: %.3e'%(dt_thermal, par.dt)
      par.dt = np.min([par.dt, dt_thermal])
   
   """
   if par.MHD:   # Maybe there is a problem with the time stepping.. dont think so
      c = np.sqrt(4.*np.pi)
      vA = np.sqrt( np.max( (var.Bx*var.Bx + var.By*var.By + var.Bz*var.Bz)  /   (var.rho*(4.*np.pi/(c*c) ) ) ) )
      v = vA/(np.sqrt(1+vA*vA/(c*c)))
      
      dt_alfven = par.cfl_set*Grid.dz/v   # We are assuming uniform grid
      #print dt_alfven
      par.dt = np.min([par.dt, dt_alfven])
   """
   par.cfl = par.dt*vcharZ/Grid.dz
