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
#   print var.rho
   c_s = np.sqrt(par.gamma*var.P[1:-1,1:-1]/var.rho[1:-1,1:-1])
   #vcharZ = np.max([np.max(abs(var.vZ + c_s)), np.max(abs(var.vZ - c_s))])
   
   if par.dim == 2 and not par.MHD:
      vcharY = np.max([np.max(abs(var.vY + c_s)), np.max(abs(var.vY - c_s))])
      par.dt = par.cfl_set/( vcharZ/Grid.dz + vcharY/Grid.dz  )
   elif par.dim==1 and not par.MHD:
      par.dt = par.cfl_set*Grid.dz/vcharZ

   if par.ThermalDiffusion and not par.ImplicitConduction:
      dt_thermal = par.f_cfl*np.min(par.cv*var.rho*Grid.dz*Grid.dz/var.kappa) 
      #print 'Thermal: %.3e \t Sound: %.3e'%(dt_thermal, par.dt)
      par.dt = np.min([par.dt, dt_thermal])
   
   if par.MHD:
      vA = np.sqrt(  (var.Bx*var.Bx + var.By*var.By + var.Bz*var.Bz)  /   var.rho )[1:-1,1:-1]   
      v = np.sqrt(var.vX*var.vX + var.vY*var.vY + var.vZ*var.vZ)[1:-1, 1:-1]
      v_ms = np.sqrt(vA**2 + c_s**2)

      vcharP = np.max( np.abs(v+v_ms) )
      vcharM = np.max( np.abs(v-v_ms) )

      vchar = np.max([vcharP,vcharM])
      if vchar!=0: 
         par.dt = par.cfl_set*Grid.dz/vchar   # We are assuming uniform grid
         #par.dt = np.min([par.dt, dt_alfven])
   
   #par.cfl = par.dt*vcharZ/Grid.dz







