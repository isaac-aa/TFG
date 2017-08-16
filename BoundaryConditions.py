#--------------------------------------------
#		BoundaryConditions.py
# This module defines the different Boundary
# conditions that can be applied to the ghost
# cells. The common argument "args" is a list
# of user defined values. Its length can depend
# on the BC
#
#--------------------------------------------


import numpy as np
import Parameters as par
import Variables as var
import Grid

print 'Loading BC..'


class BoundaryCondition(object):
   """
   Parent class for the Boundary Conditions.
   Common to all the BC is the region where they are set, so this is the same for all
   """


   def __init__(self, *args):
      """
      Initialization function.
      It can accept the region at which the BC is applied to
      """
      
      if args[0] != ():
         self.setRegion(*args[0])
      self.name = 'Generic BC'
      
   def setRegion(self, region):
      """
      The slice for the BC is done. 'region' can be:
	    'R'igth 
	    'L'eft
	    'T'op
	    'B'ottom
	    Also taken into account the number of dimensions of the problem
      """
     
      self.region = region


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


   def setStaggered(self):
      """
      If the mesh is staggered, the 'R' and 'T' boundaries are shifted
      """
      if self.region=='R':
         if par.dim == 1:
            self.sliceBC = -2
            self.sliceOne = -3
            self.sliceTwo = -4

   def __str__(self):
      return 'BC: ' + self.name + ' @ ' + self.region


   def setup(self):
      """
      Setup function to pre-allocate needed data
      """
      print 'Setup \t '+ str(self.__str__())

   def computeBC(self):
      """
      Function to be overrided by each BC
      """
      print '##WARNING## \n \n Select a BC' 








class BCComposite(BoundaryCondition):
   """
   Composite for different boundaries for each variable
   """

   name = 'BCComposite'
   
   def setup(self, rhoBC, momentumBC, energyBC, conf):
      """
         Setup for the multiple BC
         To allow for more than one dimension within the same function,
         momentumBC should be a tuple of the dimension par.dim
         conf is a tuple of size 3, storing the extra arguments for the BC if needed
      """
      
      
      print 'Setup \t ', str(self) 
      if par.dim==1:
         self.rhoBC = rhoBC(self.region)
         self.momentumZBC = momentumBC[0](self.region)
         self.momentumYBC = None
         self.energyBC = energyBC(self.region)

         if conf[0] != None: self.rhoBC.setup(conf[0])
         if conf[1] != None: self.momentumZBC.setup(conf[1])
         if conf[2] != None: self.energyBC.setup(conf[2])
         
         if par.staggered:
            self.momentumZBC.setStaggered()
                        
      elif par.dim==2:
         self.rhoBC = rhoBC(self.region)
         self.momentumZBC = momentumBC[0](self.region)
         self.momentumYBC = momentumBC[1](self.region)
         self.energyBC = energyBC(self.region)

         if conf[0] != None: self.rhoBC.setup(conf[0])
         if conf[1] != None: self.momentumZBC.setup(conf[1][0])
         if conf[2] != None: self.momentumYBC.setup(conf[1][1])
         if conf[3] != None: self.energyBC.setup(conf[2])
       
   def computeBC(self, rho, momentum, energy):
      
      self.rhoBC.computeBC(rho)
      self.momentumZBC.computeBC(momentum[0])
      if self.momentumYBC!=None: self.momentumYBC.computeBC(momentum[1])
      if par.ImplicitConduction and par.ThermalDiffusion:
         self.energyBC.computeImplicitBC()
      else:
         self.energyBC.computeBC(energy)
   



class ZeroDer(BoundaryCondition):
   name = 'Zero derivative across the boundary (Symmetry condition)'
   def computeBC(self, variable):
      variable[self.sliceBC] = variable[self.sliceOne]

class AntiSym(BoundaryCondition):
   name = 'Zero at the boundary (Antisymmetry condition)'
   def computeBC(self, variable):
         variable[self.sliceBC] = -variable[self.sliceOne]


class Fixed(BoundaryCondition):
   name = 'Fixed value'

   def setup(self, fixedVar):
      self.fixed = fixedVar

   def computeBC(self, variable):
      variable[self.sliceBC] = 2*self.fixed - variable[self.sliceOne]


class FixedT(BoundaryCondition):
   name = 'Fixed temperature value'

   def setup(self, fixedVar):
      self.fixed = fixedVar

   def computeBC(self, variable):
      variable[self.sliceBC] = 2*self.fixed - variable[self.sliceOne]
      boundaryE = 0.5*(var.rho[self.sliceBC] + var.rho[self.sliceOne])*par.cv*self.fixed 
      if par.dim == 1:
         var.energy[self.sliceBC] = 2*boundaryE - par.cv*var.rho[self.sliceOne]*var.T[self.sliceOne] + 0.5*var.momentumZ[self.sliceBC]*var.momentumZ[self.sliceBC]/var.rho[self.sliceBC]
      elif par.dim == 2:
         var.energy[self.sliceBC] = 2*boundaryE - par.cv*var.rho[self.sliceOne]*var.T[self.sliceOne] + 0.5*(var.momentumZ[self.sliceBC]*var.momentumZ[self.sliceBC]+ var.momentumY[self.sliceBC]*var.momentumY[self.sliceBC])/var.rho[self.sliceBC]

   def computeImplicitBC(self):
      self.computeBC(var.energy)
      var.diag[self.sliceBC] = 1.    #This only works for one-dimensional cases
      if self.region=="L":
         var.upper[1] = 1.
      if self.region=="R":
         var.lower[-2] = 1.
      var.rhs[self.sliceBC] = 2.*var.rho[self.sliceBC]*par.cv*self.fixed


class ConstantSecondDer(BoundaryCondition):
   name = 'Fixed second derivative'

   def computeBC(self, variable):
      variable[self.sliceBC] = 2*variable[self.sliceOne] - variable[self.sliceTwo]



class Hydrostatic(BoundaryCondition):
   name = 'Hydrostatic condition'

   def computeBC(self, variable):
      if par.dim == 1:
         dlogP = Grid.dz*0.5*(var.rho[self.sliceBC] + var.rho[self.sliceOne])*np.abs(par.g)
         if self.region == "R":
            var.P[self.sliceBC] = np.exp( np.log(var.P[self.sliceOne]) + dlogP )
         if self.region == "L":
            var.P[self.sliceBC] = np.exp( np.log(var.P[self.sliceOne]) - dlogP ) 
         var.energy[self.sliceBC] = var.P[self.sliceBC]/(par.gamma-1.) + 0.5*var.momentumZ[self.sliceBC]*var.momentumZ[self.sliceBC]/var.rho[self.sliceBC]
      elif par.dim == 2:
               var.energy[self.sliceBC] = var.P[self.sliceBC]/(par.gamma-1.) + 0.5*(var.momentumZ[self.sliceBC]*var.momentumZ[self.sliceBC]+ var.momentumY[self.sliceBC]*var.momentumY[self.sliceBC])/var.rho[self.sliceBC]


   def computeImplicitBC(self):
      self.computeBC(var.energy)
      PJump = Grid.dz*0.5*(var.rho[self.sliceBC]+var.rho[self.sliceOne])*np.abs(par.g)
      if self.region=="L":
         var.diag[self.sliceBC] = 1.
         var.upper[1] = -1.
      if self.region=="R":      
         var.diag[self.sliceBC] = 1.
         var.lower[-2] = -1.

      var.rhs[self.sliceBC] = PJump/(par.gamma-1.)




class Periodic(BoundaryCondition):
   """
    Given two generic BC region, it connects them
   """
   
   name = 'Periodic'

   def __init__(self, regionA, regionB):
      self.regionA = regionA
      self.regionB = regionB

   def __str__(self):
      return 'Periodic BC: ' +  self.regionA.region + ' & ' + self.regionB.region 

   def computeBC(self, rho, momentum, energy):
      rho[self.regionA.sliceBC] = rho[self.regionB.sliceOne]
      rho[self.regionB.sliceBC] = rho[self.regionA.sliceOne]
      
      momentum[0][self.regionA.sliceBC] = momentum[0][self.regionB.sliceOne]
      momentum[0][self.regionB.sliceBC] = momentum[0][self.regionA.sliceOne]

      if par.dim == 2:
         momentum[1][self.regionA.sliceBC] = momentum[1][self.regionB.sliceOne]
         momentum[1][self.regionB.sliceBC] = momentum[1][self.regionA.sliceOne]

      energy[self.regionA.sliceBC] = energy[self.regionB.sliceOne]
      energy[self.regionB.sliceBC] = energy[self.regionA.sliceOne]







