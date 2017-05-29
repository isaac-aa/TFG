#--------------------------------------------
#		Main.py
# This module is the core of the simulation.
# Firstly it imports all the neccesary modules,
# creates the mesh, set the initial state and
# starts the main loop
#
#--------------------------------------------

import numpy as np
import time

CPUTime0 = time.clock()

import Parameters as par
import Grid

# ------------------ MESH CREATION -------------
#Grid.Uniform1DGrid(par.N, par.z0, par.zf)
Grid.ReadGridFromFile('Extras/ThermalEq_IC_2.dat') #('hydrostatic_equilibrium_2.dat')

import Variables as var
import Settings as sets
import SourceTerm
import TimeStep
import ChangeOfVar
import Save

# ------------------ INITIAL CONDITION ---------

sets.InitialCondition(sets.argsIC)
ChangeOfVar.ConvertToPrim()
Save.Plot()
# ------------------ MAIN LOOP -----------------

CPUTimeConf = time.clock()
print '################ \n'
print 'Configuration took %.4f s'%(CPUTimeConf-CPUTime0)
print '################ \t SIMULATION START \n'

WallTime0 = time.time()

while (par.it<=par.max_it and par.tt<=par.tf):
   TimeStep.ComputeDT()  

   par.tt += par.dt  
   par.it += 1
   var.lastrho[:] = var.rho[:]
   var.lastmomentum[:] = var.momentum[:]
   var.lastenergy[:] = var.energy[:]

   # Time step
   sets.Scheme() 
   
   # Source computation
   if par.IsComputingSource:
      SourceTerm.ComputeSource()

   # Boundary conditions     
   sets.BoundaryConditionL(sets.argsL)
   sets.BoundaryConditionR(sets.argsR)

   # Compute change of variables
   ChangeOfVar.ConvertToPrim()

   if par.it%par.save_rate == 0.:
     ItTime = time.clock()
     print ' ## IT: %d \t CPU-Time: %.2f s \t Wall-time: %.2f'%(par.it, ItTime-CPUTimeConf, time.time()-WallTime0)
     print 'DT: %.3e \t t: %.3f s \t CFL: %.3e'%(par.dt, par.tt, par.cfl)
     Save.Plot()
     maxIndex = np.argmax(np.abs(var.v))
     print 'Max v: %.2f cm/s @ z[%d]: %.3e \n'%(var.v[maxIndex], maxIndex, Grid.z[maxIndex])
     #print var.P[0], var.P[1], var.P[2], 0.5*Grid.dz*boundaryRho*np.abs(par.g)







