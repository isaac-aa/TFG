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
import Grid

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
   
   if par.ImplicitConduction:
      var.diag[i] = 1.    
      if args[0]=="L":
         var.upper[1] = 1.
      if args[0]=="R":
         var.lower[-2] = 1. 
      var.rhs[i] = 2.*var.rho[i_one]*par.cv*var.T[i_one]      
   else:
      var.energy[i] = var.energy[i_one]


def Wall2D(args):
   if args[0] == 'R':
      var.rho[:,-1] = var.rho[:,-2]
      var.momentumZ[:,-1] = var.momentumZ[:,-2]
      var.momentumY[:,-1] = var.momentumY[:,-2]
      var.energy[:,-1] = var.energy[:,-2]


def WallFixedRhoHydrostaticP2D(args):
   rhoFixed = args[1]
   # USE MASKED ARRAYS
   if args[0]=='L':
      var.rho[:, 0] = 2.*rhoFixed - var.rho[:,1]  
      var.P[:, 0] = var.P[:,1] + 0.5*Grid.dz*rhoFixed*np.abs(par.g)
      var.momentumZ[:,0] = var.momentumZ[:,1]
      var.momentumY[:,0] = var.momentumY[:,1]
      var.energy[:,0] = var.P[:,0]/(par.gamma-1.) + 0.5*(var.momentumZ[:,0]*var.momentumZ[:,0] - var.momentumY[:,0]*var.momentumY[:,0])/var.rho[:,0]

   if args[0]=='R':
      var.rho[:,-1] = 2.*rhoFixed - var.rho[:,-2]
      var.P[:, -1] = var.P[i_one] + 0.5*Grid.dz*rhoFixed*np.abs(par.g)
   if args[0]=='T':
      var.rho[-1,:] = 2.*rhoFixed - var.rho[-2,:]
      var.P[-1,:] = var.P[i_one] + 0.5*Grid.dy*rhoFixed*np.abs(par.g)
   if args[0]=='B':
      var.rho[0,:] = 2.*rhoFixed - var.rho[1,:]
      var.P[0,:] = var.P[i_one] + 0.5*Grid.dy*rhoFixes*np.abs(par.g)


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
   
   FixedT = args[1]
   var.rho[i] = var.rho[i_one]
   
   var.momentum[i] = -var.momentum[i_one]
   
   #E_k = 0.5*var.v[i_one]*var.v[i_one]*var.rho[i_one]
   #boundaryE = var.rho[i]*par.cv*args[1] + E_k
   #var.energy[i] = 2*boundaryE-var.energy[i_one]
   if par.ImplicitConduction:
      var.diag[i] = 1.    
      if args[0]=="L":
         var.upper[1] = 1.
      if args[0]=="R":
         var.lower[-2] = 1. 
      var.rhs[i] = 2.*var.rho[i]*par.cv*FixedT
   else:
      boundaryE = var.rho[i_one]*par.cv*FixedT 
      internalE = 2*boundaryE - par.cv*var.rho[i_one]*var.T[i_one]   
      var.energy[i] = internalE + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]

def WallSecondRhoFixedT(args):
   if args[0]=="L":
     i = 0
     i_one = 1
     i_two = 2
   if args[0]=="R":
     i = -1
     i_one = -2
     i_two = -3
   
   FixedT = args[1]
   
   var.rho[i] = 2.*var.rho[i_one]-var.rho[i_two]
   boundaryRho = 0.5*(var.rho[i]+var.rho[i_one])
   
   var.momentum[i] = -var.momentum[i_one]    #v = 0
   
   if par.ImplicitConduction:
      var.diag[i] = 1.    
      if args[0]=="L":
         var.upper[1] = 1.
      if args[0]=="R":
         var.lower[-2] = 1. 
      var.rhs[i] = 2.*boundaryRho*par.cv*FixedT
   else:
      boundaryE = boundaryRho*par.cv*FixedT 
      internalE = 2*boundaryE - par.cv*var.rho[i_one]*var.T[i_one]   
      var.energy[i] = internalE + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]


def SymVSecondRhoFixedT(args):
   if args[0]=="L":
     i = 0
     i_one = 1
     i_two = 2
   if args[0]=="R":
     i = -1
     i_one = -2
     i_two = -3

   FixedT = args[1]

   var.rho[i] = 2.*var.rho[i_one]-var.rho[i_two]
   boundaryRho = 0.5*(var.rho[i]+var.rho[i_one])

   var.momentum[i] = var.momentum[i_one]    #v = 0

   if par.ImplicitConduction:
      var.diag[i] = 1.
      if args[0]=="L":
         var.upper[1] = 1.
      if args[0]=="R":
         var.lower[-2] = 1.
      var.rhs[i] = 2.*boundaryRho*par.cv*FixedT
   else:
      boundaryE = boundaryRho*par.cv*FixedT
      internalE = 2*boundaryE - par.cv*var.rho[i_one]*var.T[i_one]
      var.energy[i] = internalE + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]



def SymVFixedRhoHydrostaticP(args):
   if args[0]=="L":
     i = 0
     i_one = 1
     i_two = 2
   if args[0]=="R":
     i = -1
     i_one = -2
     i_two = -3

   FixedRho = args[1]

   var.rho[i] = 2.*FixedRho-var.rho[i_one]
   boundaryRho = FixedRho

   var.momentum[i] = var.momentum[i_one]   

   if par.ImplicitConduction:
      PJump = Grid.dz*boundaryRho*np.abs(par.g)
      if args[0]=="L":
         var.diag[i] = 1.
         var.upper[1] = -1.
      if args[0]=="R":      #This case is not used
         var.diag[i] = 1.
         var.lower[-2] = -1.

      var.rhs[i] = PJump/(par.gamma-1.)
   else:
      var.P[i] = var.P[i_one] + 0.5*Grid.dz*boundaryRho*np.abs(par.g)
      var.energy[i] = var.P[i]/(par.gamma-1.) + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]



def WallFixedRhoHydrostaticP(args):
   if args[0]=="L":
     i = 0
     i_one = 1
     i_two = 2
   if args[0]=="R":
     i = -1
     i_one = -2
     i_two = -3

   
   FixedRho = args[1]

   var.rho[i] = 2.*FixedRho-var.rho[i_one]
   boundaryRho = 0.5*(var.rho[i]+var.rho[i_one])

   var.momentum[i] = -var.momentum[i_one]    #v = 0

   if par.ImplicitConduction:
      PJump = Grid.dz*boundaryRho*np.abs(par.g)
      if args[0]=="L":
         var.diag[i] = 1.
         var.upper[1] = -1.
      if args[0]=="R":      #This case is not used
         var.diag[i] = 1.
         var.lower[-2] = -1.

      var.rhs[i] = PJump/(par.gamma-1.)
   else:
      dlogP = Grid.dz*boundaryRho*np.abs(par.g)
      if args[0] == "R":
          var.P[i] = np.exp( np.log(var.P[i_one]) + dlogP )
      if args[1] == "L":
          var.P[i] = np.exp( np.log(var.P[i_one]) - dlogP ) 
      var.energy[i] = var.P[i]/(par.gamma-1.) + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]





def WallSecondRhoHydrostaticP(args):
   if args[0]=="L":
     i = 0
     i_one = 1
     i_two = 2
   if args[0]=="R":
     i = -1
     i_one = -2
     i_two = -3

   var.rho[i] = 2.*var.rho[i_one]-var.rho[i_two]
   boundaryRho = 0.5*(var.rho[i]+var.rho[i_one])

   var.momentum[i] = -var.momentum[i_one]    #v = 0

   if par.ImplicitConduction:
      PJump = Grid.dz*boundaryRho*np.abs(par.g)
      if args[0]=="L":
         var.diag[i] = 1. 
         var.upper[1] = -1.
      if args[0]=="R":      #This case is not used
         var.diag[i] = 1. 
         var.lower[-2] = -1.

      var.rhs[i] = PJump/(par.gamma-1.)
   else:
      dlogP = Grid.dz*boundaryRho*np.abs(par.g)
      if args[0] == "L":
          var.P[i] = np.exp( np.log(var.P[i_one]) + dlogP )
      if args[1] == "R":
          var.P[i] = np.exp( np.log(var.P[i_one]) - dlogP )

      var.energy[i] = var.P[i]/(par.gamma-1.) + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]


def SymVSecondRhoHydrostaticP(args):
   if args[0]=="L":
     i = 0
     i_one = 1
     i_two = 2
   if args[0]=="R":
     i = -1
     i_one = -2
     i_two = -3

   var.rho[i] = 2.*var.rho[i_one]-var.rho[i_two]
   boundaryRho = 0.5*(var.rho[i]+var.rho[i_one])

   var.momentum[i] = var.momentum[i_one]    #Symmetric momentum

   if par.ImplicitConduction:
      PJump = Grid.dz*boundaryRho*np.abs(par.g)
      if args[0]=="L":
         var.diag[i] = 1.
         var.upper[1] = -1.
      if args[0]=="R":      #This case is not used
         var.diag[i] = 1.
         var.lower[-2] = -1.

      var.rhs[i] = PJump/(par.gamma-1.)
   else:
      var.P[i] = var.P[i_one] + 0.5*Grid.dz*boundaryRho*np.abs(par.g)
      var.energy[i] = var.P[i]/(par.gamma-1.) + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]


def WallFixedRhoFixedT(args):
   if args[0]=="L":
     i = 0
     i_one = 1
   if args[0]=="R":
     i = -1
     i_one = -2
   
   FixedT = args[1]
   FixedRho = args[2]
   
   var.rho[i] = FixedRho #2.*FixedRho-var.rho[i_one]
   
   var.momentum[i] = -var.momentum[i_one]    #v = 0

      
   if par.ImplicitConduction:
      var.diag[i] = 1.    
      if args[0]=="L":
         var.upper[1] = 0. #1.
      if args[0]=="R":
         var.lower[-2] = 0. #1. 
      var.rhs[i] = FixedRho*par.cv*FixedT #2.*FixedRho*par.cv*FixedT
   else:
      boundaryE = FixedRho*par.cv*FixedT 
      internalE = 2*boundaryE - par.cv*var.rho[i_one]*var.T[i_one]   
      var.energy[i] = internalE + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]   


def SymVFixedRhoFixedT(args):
   if args[0]=="L":
     i = 0
     i_one = 1
   if args[0]=="R":
     i = -1
     i_one = -2

   FixedT = args[1]
   FixedRho = args[2]

   var.rho[i] = FixedRho #2.*FixedRho-var.rho[i_one]

   var.momentum[i] = var.momentum[i_one]    #v = 0


   if par.ImplicitConduction:
      var.diag[i] = 1.
      if args[0]=="L":
         var.upper[1] = 0. #1.
      if args[0]=="R":
         var.lower[-2] = 0. #1. 
      var.rhs[i] = FixedRho*par.cv*FixedT #2.*FixedRho*par.cv*FixedT
   else:
      boundaryE = FixedRho*par.cv*FixedT
      internalE = 2*boundaryE - par.cv*var.rho[i_one]*var.T[i_one]
      var.energy[i] = internalE + 0.5*var.momentum[i]*var.momentum[i]/var.rho[i]



   
def Periodic(args):
   var.rho[0] = var.rho[-2]
   var.rho[-1] = var.rho[1]

   var.momentum[0] = var.momentum[-2]
   var.momentum[-1] = var.momentum[1]

   var.energy[0] = var.energy[-2]
   var.energy[-1] = var.energy[1]


def Periodic2D(args):
  var.rho[:, 0] = var.rho[:,-2]
  var.rho[:, -1] = var.rho[:,1]
  var.rho[0,:] = var.rho[-2,:]
  var.rho[-1,:] = var.rho[1,:]

  
  var.momentumZ[:, 0] = var.momentumZ[:,-2]
  var.momentumZ[:, -1] = var.momentumZ[:,1]
  var.momentumZ[0,:] = var.momentumZ[-2,:]
  var.momentumZ[-1,:] = var.momentumZ[1,:]
  
  var.momentumY[:, 0] = var.momentumY[:,-2]
  var.momentumY[:, -1] = var.momentumY[:,1]
  var.momentumY[0,:] = var.momentumY[-2,:]
  var.momentumY[-1,:] = var.momentumY[1,:]
 
  var.energy[:, 0] = var.energy[:,-2]
  var.energy[:, -1] = var.energy[:,1]
  var.energy[0,:] = var.energy[-2,:]
  var.energy[-1,:] = var.energy[1,:]



def PeriodicY2D(args):
  var.rho[0,:] = var.rho[-2,:]
  var.rho[-1,:] = var.rho[1,:]

  var.momentumZ[0,:] = var.momentumZ[-2,:]
  var.momentumZ[-1,:] = var.momentumZ[1,:]
  
  var.momentumY[0,:] = var.momentumY[-2,:]
  var.momentumY[-1,:] = var.momentumY[1,:]
 
  var.energy[0,:] = var.energy[-2,:]
  var.energy[-1,:] = var.energy[1,:]



