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


class BoundaryCondition(object):
   """
   Parent class for the Boundary Conditions
   """


   def __init__(self, region):
     """
     The mask for the BC is done. Options:
	'R'igth 
	'L'eft
	'T'op
	'B'ottom
     """
     self.region = region
     self.name = 'Generic BC'


     if region=='R':
        if par.dim == 2:
          self.sliceBC = ( slice(None,None), -1)
          self.sliceOne = ( slice(None,None), -2)
          self.sliceTwo = ( slice(None,None), -3)
        if par.dim == 1:
          self.sliceBC = -1
          self.sliceOne = -2
          self.sliceTwo = -3
     if region=='L':
        if par.dim == 2:
          self.sliceBC = ( slice(None,None), 0)
          self.sliceOne = ( slice(None,None), 1)
          self.sliceTwo = ( slice(None,None), 2)
        if par.dim == 1:
          self.sliceBC = 0
          self.sliceOne = 1
          self.sliceTwo = 2
     if region=='T':
        if par.dim == 2:
          self.sliceBC = ( -1, slice(None,None))
          self.sliceOne = ( -2, slice(None,None))
          self.sliceTwo = ( -3, slice(None,None))
     if region=='B':
        if par.dim == 2:	
          self.sliceBC = ( 0, slice(None,None))
          self.sliceOne = ( 1, slice(None,None))
          self.sliceTwo = ( 2, slice(None,None))


   def __str__(self):
     return 'BC: ' + self.name + ' @ ' + self.region


   def setup(self):
     """
     Setup function to pre-allocate neccesary data
     """
     print 'Setup \t '+ str(self.__str__())

   def computeBC(self):
     """
     Function to be overrided by each BC
     """
     print 'Select a BC' 








class BCComposite(BoundaryCondition):
   """
   Composite for different boundaries for each variable
   """

   name = 'BCComposite'
   
   def setup(self, rhoBC, momentumBC, energyBC):
      """ Setup for a one-dimensional BC """
      print 'Setup \t ', self 
      self.rhoBC = rhoBC
      self.momentumZBC = momentumBC
      self.momentumYBC = None
      self.energyBC = energyBC

      self.rhoBC.setup()
      self.momentumZBC.setup()
      self.energyBC.setup()

   def setup(self, rhoBC, momentumZBC, momentumYBC, energyBC):
      """ Setup for a two-dimensional BC """
      self.rhoBC = rhoBC
      self.momentumZBC = momentumZBC
      self.momentumYBC = momentumYBC
      self.energyBC = energyBC

      self.rhoBC.setup()
      self.momentumZBC.setup()
      self.momentumYBC.setup()
      self.energyBC.setup()


   def computeBC(self):
      rhoBC.computeBC(var.rho)
      momentumZBC.computeBC(var.momentum)
      if momentumYBC!=None: momentumYBC.computeBC(var.momentum)
      if par.ImplicitCondution:
        energyBC.computeImplicitBC(var.energy)
      else:
        energyBC.computeBC(var.energy)
   



class ZeroDer(BoundaryCondition):
   name = 'Zero der'
   def computeBC(self, variable):
     variable[self.sliceBC] = variable[self.sliceOne]

class AntiSym(BoundaryCondition):
   name = 'Antisymmetry'
   def computeBC(self, variable):
     variable[self.sliceBC] = -variable[self.sliceOne]


class Fixed(BoundaryCondition):
   name = 'Fixed'

   def setup(self, fixedVar):
     self.fixed = fixedVar


   def computeBC(self, variable):
     variable[sliceBC] = 2*self.fixed - variable[sliceOne]

   def computeImplicitBC(self):
      var.diag[i] = 1.
      if args[0]=="L":
         var.upper[1] = 1.
      if args[0]=="R":
         var.lower[-2] = 1.
      var.rhs[i] = 2.*var.rho[i]*par.cv*self.fixed


class ConstantSecondDer(BoundaryCondition):
   name = 'Fixed second derivative'

   def computeBC(self, variable):
      variable[sliceBC] = 2*variable[sliceOne] - variable[sliceTwo]



class Hydrostatic(BoundaryCondition):
   name = 'Hydrostatic condition'

   def computeBC(self, variable):
      var.energy[sliceBC] = var.P[sliceBC]/(par.gamma-1.) + 0.5*var.momentum[sliceBC]*var.momentum[sliceBC]/var.rho[sliceBC]


   def computeImplicitBC(self):
      PJump = Grid.dz*boundaryRho*np.abs(par.g)
      if args[0]=="L":
         var.diag[i] = 1.
         var.upper[1] = -1.
      if args[0]=="R":      
         var.diag[i] = 1.
         var.lower[-2] = -1.

      var.rhs[i] = PJump/(par.gamma-1.)




class Periodic(BoundaryCondition):
   """ Given two generic BC region, it connects them """
   name = 'Periodic'

   def __init__(self, regionA, regionB):
     self.regionA = regionA
     self.regionB = regionB

   def __str__(self):
     return 'Periodic BC: ' +  str(self.regionA) + str(self.regionB) 

   def computeBC(self):
     var.rho[self.regionA.sliceBC] = var.rho[self.regionB.sliceOne]
     var.rho[self.regionB.sliceBC] = var.rho[self.regionA.sliceOne]

     var.momentum[self.regionA.sliceBC] = var.momentum[self.regionB.sliceOne]
     var.momentum[self.regionB.sliceBC] = var.momentum[self.regionA.sliceOne]

     var.energy[self.regionA.sliceBC] = var.energy[self.regionB.sliceOne]
     var.energy[self.regionB.sliceBC] = var.energy[self.regionA.sliceOne]



# ---------------- OLD PART
"""
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


"""


