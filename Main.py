import numpy as np
import time


import Variables as var
import Settings as sets
import SourceTerm
import TimeStep
import ChangeOfVar
import Save

# ------------------ TIME STEP -----------------

TimeStep.ComputeDT()

while (var.it<=var.max_it and var.tt<=var.tf):

   # Boundary conditions
   sets.BoundaryConditionL()
   sets.BoundaryConditionR()

   # Source computation
   if var.IsComputingSource:
      SourceTerm.ComputeSource()

   # Time step
   sets.Scheme() 

   # Compute change of variables
   ChangeOfVar.ConvertToPrim()

   TimeStep.ComputeDT()


   if var.it%var.save_rate == 0.:
     Save.Plot()
     print var.it, var.tt, var.dt, var.cfl

   var.tt += var.dt
   var.it += 1
   #time.sleep(1)






