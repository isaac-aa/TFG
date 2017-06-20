#--------------------------------------------
#		TimeStep.py
# This module in on charge of computing the 
# time step each iteration
#
#--------------------------------------------

import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
comm_size = comm.Get_size()

import Parameters as par
import Grid
import Variables as var

print "Loading TimeStep.."

def ComputeDT():
   c_s = np.sqrt(par.gamma*var.P/var.rho)
   vchar1 = np.max(abs(var.v + c_s))
   vchar2 = np.max(abs(var.v - c_s))
   
   par.dt = par.cfl_set*Grid.dz/np.max( [vchar1, vchar2] )

   if par.ThermalDiffusion and not par.ImplicitConduction:
      dt_thermal = par.f_cfl*np.min(par.cv*var.rho*Grid.dz*Grid.dz/var.kappa) 
      #print 'Thermal: %.3e \t Sound: %.3e'%(dt_thermal, par.dt)
      par.dt = np.min([par.dt, dt_thermal])
      
   par.cfl = par.dt*np.max( [vchar1, vchar2] )/Grid.dz

   comm.Barrier()
   all_dt = None
   if rank==0:
     all_dt = np.zeros([comm_size])
   comm.Gather(par.dt,all_dt, root=0)
   if rank==0:
     min_dt = np.min(all_dt)
   else:
     min_dt = None
   min_dt = comm.bcast(min_dt, root=0)
   par.dt = min_dt
