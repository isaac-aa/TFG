import numpy as np
import Variables as var
import Settings as sets

print 'Loading BC..'

def Wall():
   if sets.args[0]=="L":
     i = 0
     i_one = 1
   if sets.args[0]=="R":
     i = -1
     i_one = -2

   var.rho[i] = var.rho[i_one]
   var.momentum[i] = var.momentum[i_one]
   var.energy[i] = var.energy[i_one]

def FixedRhoP():
   if sets.args[0]=="L":
     i = 0
     i_one = 1
   if sets.args[0]=="R":
     i = -1
     i_one = -2

   var.rho[i] = sets.args[1] 
   var.momentum[i] = var.momentum[i_one]
   var.energy[i] = sets.args[2]/(sets.args[1]*(var.gamma-1.) ) #energy[i_one]

def Periodic():
   var.rho[0] = var.rho[-2]
   var.rho[-1] = var.rho[1]

   var.momentum[0] = var.momentum[-2]
   var.momentum[-1] = var.momentum[1]

   var.energy[0] = var.energy[-2]
   var.energy[-1] = var.energy[1]





