import numpy as np

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
   energy[i] = energy[i_one]

"""
def SoundWave(side, rho, momentum, energy, args):  #args=[t]
   rho[0] = rho[1]
   rho[-1] = rho[-2]
   momentum[-1] = 0.2*np.sin(2*np.pi*t/0.1)
   momentum[0] = 0.
   energy[0] = energy[1]
   energy[-1] = energy[-2]
"""
