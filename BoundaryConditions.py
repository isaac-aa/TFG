import numpy as np
from Variables import *
print 'Loading BC..'

def Wall(rho, momentum, energy, args):
   if args[0]=="L":
     i = 0
     i_one = 1
   if args[0]=="R":
     i = -1
     i_one = -2

   rho[i] = rho[i_one]
   momentum[i] = momentum[i_one]
   energy[i] = energy[i_one]

def FixedRho(rho, momentum,energy,args):
   if args[0]=="L":
     i = 0
     i_one = 1
   if args[0]=="R":
     i = -1
     i_one = -2

   rho[i] = args[1] 
   momentum[i] = momentum[i_one]
   energy[i] = 1./(gamma-1.) #energy[i_one]

def Periodic(rho, momentum, energy, args):
   rho[0] = rho[-2]
   rho[-1] = rho[1]

   momentum[0] = momentum[-2]
   momentum[-1] = momentum[1]

   energy[0] = energy[-2]
   energy[-1] = energy[1]





