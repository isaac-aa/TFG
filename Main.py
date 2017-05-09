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


import Parameters as par
import Grid

# ------------------ MESH CREATION -------------
Grid.Uniform1DGrid(par.N, par.z0, par.zf)

import Variables as var
import Settings as sets
import SourceTerm
import TimeStep
import ChangeOfVar
import Save

# ------------------ INITIAL CONDITION ---------

sets.InitialCondition(sets.argsIC)
ChangeOfVar.ConvertToPrim()

# ------------------ MAIN LOOP -----------------

   
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
   # NOTE: This should be done after the Scheme if doing Operator Splitting (TODO)
   if par.IsComputingSource:
      SourceTerm.ComputeSource()

   # Boundary conditions     
   sets.BoundaryConditionL(sets.argsL)
   sets.BoundaryConditionR(sets.argsR)

   # Compute change of variables
   ChangeOfVar.ConvertToPrim()

   if par.it%par.save_rate == 0.:
     Save.Plot()
     print par.it, par.tt, par.dt, par.cfl








