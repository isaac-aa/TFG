#--------------------------------------------
#		BoundaryConditions.py
# This module defines the different Boundary
# conditions that can be applied to the ghost
# cells. The common argument "args" is a list
# of user defined values. Its lenght can depend
# on the BC
#
#--------------------------------------------


import numpy as np
import Parameters as par
import Variables as var
import Settings as sets

print 'Loading BC..'

def Wall(args):
   if args[0]=="L":
     i = 0
     i_one = 1
   if args[0]=="R":
     i = -1
     i_one = -2

   var.rho[i] = var.rho[i_one]
   var.momentum[i] = var.momentum[i_one]
   var.energy[i] = var.energy[i_one]

def FixedRhoP(args):
   if args[0]=="L":
     i = 0
     i_one = 1
   if args[0]=="R":
     i = -1
     i_one = -2

   var.rho[i] = args[1] 
   var.momentum[i] = var.momentum[i_one]
   var.energy[i] = args[2]/(par.gamma-1.) + 0.5*var.rho[i]*var.v[i_one]*var.v[i_one] #energy[i_one]
   
def FixedT(args):
   if args[0]=="L":
     i = 0
     i_one = 1
   if args[0]=="R":
     i = -1
     i_one = -2
   
   var.rho[i] = var.rho[i_one]
   var.momentum[i] = var.momentum[i_one]
   E_k = 0.5*var.v[i_one]*var.v[i_one]*var.rho[i]
   var.energy[i] = var.rho[i]*par.cv*args[1] + E_k

def Periodic(args):
   var.rho[0] = var.rho[-2]
   var.rho[-1] = var.rho[1]

   var.momentum[0] = var.momentum[-2]
   var.momentum[-1] = var.momentum[1]

   var.energy[0] = var.energy[-2]
   var.energy[-1] = var.energy[1]





