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


   def setStaggered(self, momentum):
      """
      If the mesh is staggered, the 'R' and 'T' boundaries are shifted
      """
      
      if momentum=='Z':
         if self.region=='R':
            if par.dim==1:
               self.sliceBC = -2
               self.sliceOne = -3
               self.sliceTwo = -4
            else:
               self.sliceBC = ( slice(None,None), -2)
               self.sliceOne = ( slice(None,None), -3)
               self.sliceTwo = ( slice(None,None), -4)       
         elif self.region=='L':
            if par.dim==2:
               self.sliceBC = ( slice(None,None), 0)
               self.sliceOne = ( slice(None,None), 1)
               self.sliceTwo = ( slice(None,None), 2)    
         elif self.region=='T':
            self.sliceBC = ( -2, slice(None,None))
            self.sliceOne = ( -3, slice(None,None))
            self.sliceTwo = ( -4, slice(None,None))     
         elif self.region=='B':
            self.sliceBC = ( 1, slice(None,None))
            self.sliceOne = ( 2, slice(None,None))
            self.sliceTwo = ( 3, slice(None,None))    
      
      if momentum=='Y':
         if self.region=='R':
            self.sliceBC = ( slice(None,None), -1)
            self.sliceOne = ( slice(None,None), -2)
            self.sliceTwo = ( slice(None,None), -3)              
         elif self.region=='L':
            self.sliceBC = ( slice(None,None), 1)
            self.sliceOne = ( slice(None,None), 2)
            self.sliceTwo = ( slice(None,None), 3)
         elif self.region=='T':
            self.sliceBC = ( -2, slice(None,None))
            self.sliceOne = ( -3, slice(None,None))
            self.sliceTwo = ( -4, slice(None,None))     
         elif self.region=='B':
            self.sliceBC = ( 0, slice(None,None))
            self.sliceOne = ( 1, slice(None,None))
            self.sliceTwo = ( 2, slice(None,None))                                     
                                          
      """
      if self.region=='R':
         if par.dim == 1:
            self.sliceBC = -2
            self.sliceOne = -3
            self.sliceTwo = -4
         else:
            self.sliceBC = ( slice(None,None), -2)
            self.sliceOne = ( slice(None,None), -3)
            self.sliceTwo = ( slice(None,None), -4)           
      if self.region=='T':
         if par.dim == 1:
            self.sliceBC = -2
            self.sliceOne = -3
            self.sliceTwo = -4
         else:
            self.sliceBC = ( -2, slice(None,None))
            self.sliceOne = ( -3, slice(None,None))
            self.sliceTwo = ( -4, slice(None,None))          

      """
   
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
      
   def computeImplicitBC(self):
      print '##WARNING## \n \n ' + str(self) + ' does not have implicit BC implemented'








class BCComposite(BoundaryCondition):
   """
   Composite for different boundaries for each variable
   """

   name = 'BCComposite'
   
   def setup(self, rhoBC, momentumBC, energyBC, magBC, conf):
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
         if conf[1] != None: self.momentumZBC.setup(conf[1])
         if conf[2] != None: self.momentumYBC.setup(conf[2])
         if conf[3] != None: self.energyBC.setup(conf[3])
         
         if par.staggered:
            self.momentumZBC.setStaggered('Z')
            self.momentumYBC.setStaggered('Y')
            
         if par.MHD:
            self.momentumXBC = momentumBC[2](self.region)
            self.BzBC = magBC[0](self.region)
            self.ByBC = magBC[1](self.region)
            self.BxBC = magBC[2](self.region)
      
   
   def computeBC(self, rho, momentum, energy, B):
      
      self.rhoBC.computeBC(rho)
      self.momentumZBC.computeBC(momentum[0])
      if self.momentumYBC!=None: self.momentumYBC.computeBC(momentum[1])
      if par.ImplicitConduction and par.ThermalDiffusion:
         self.energyBC.computeImplicitBC()
      else:
         self.energyBC.computeBC(energy)
   
      if par.MHD:
         if self.momentumXBC!=None: self.momentumXBC.computeBC(momentum[2])
         self.BzBC.computeBC(B[0])
         self.ByBC.computeBC(B[1])
         self.BxBC.computeBC(B[2])



class ZeroDer(BoundaryCondition):
   name = 'Zero derivative across the boundary (Symmetry condition)'
   def computeBC(self, variable):
      variable[self.sliceBC] = variable[self.sliceOne]
      
   def computeImplicitBC(self):
      self.computeBC(var.energy)
      if par.dim==2:
         rhs = np.zeros_like(var.energy[self.sliceBC])
         nnzC = 1 #- par.cv*var.rho[self.sliceOne]*var.T[self.sliceOne]
         nnzBC = -1
         return rhs, nnzC, nnzBC


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
      internalE = var.rho[self.sliceBC]*par.cv*var.T[self.sliceBC]
      boundaryE = 0.5*(var.rho[self.sliceBC] + var.rho[self.sliceOne])*par.cv*self.fixed
      
      #Not completely accurate.. should do the means for the momentums if par.staggered
      if par.staggered:
         if self.region=='R':
            Ek = 0.5*(var.momentumZ[:,-2]*var.momentumZ[:,-2] + var.momentumY[self.sliceBC]*var.momentumY[self.sliceBC])/var.rho[self.sliceBC]
         elif self.region=='T':
            Ek = 0.5*(var.momentumZ[self.sliceBC]*var.momentumZ[self.sliceBC] + var.momentumY[:,-2]*var.momentumY[:,-2])/var.rho[self.sliceBC]
         else:   
            Ek = 0.5*(var.momentumZ[self.sliceBC]*var.momentumZ[self.sliceBC] + var.momentumY[self.sliceBC]*var.momentumY[self.sliceBC])/var.rho[self.sliceBC]
      else:
         Ek = 0.5*(var.momentumZ[self.sliceBC]*var.momentumZ[self.sliceBC] + var.momentumY[self.sliceBC]*var.momentumY[self.sliceBC])/var.rho[self.sliceBC]
         
      if par.dim == 1:
         var.energy[self.sliceBC] = 2*boundaryE - par.cv*var.rho[self.sliceOne]*var.T[self.sliceOne] + Ek
      elif par.dim == 2:
         var.energy[self.sliceBC] = 2*boundaryE - par.cv*var.rho[self.sliceOne]*var.T[self.sliceOne] + Ek


   def computeImplicitBC(self):
      if par.dim == 1:
         self.computeBC(var.energy)
         var.diag[self.sliceBC] = 1.    #This only works for one-dimensional cases
         if self.region=="L":
            var.upper[1] = 1.
         if self.region=="R":
            var.lower[-2] = 1.
         var.rhs[self.sliceBC] = 2.*var.rho[self.sliceBC]*par.cv*self.fixed
      
      elif par.dim==2:
         boundaryE = 0.5*(var.rho[self.sliceBC] + var.rho[self.sliceOne])*par.cv*self.fixed
         rhs = boundaryE 
         nnzC = .5 #- par.cv*var.rho[self.sliceOne]*var.T[self.sliceOne]
         nnzBC = .5
         return rhs, nnzC, nnzBC


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

   def computeBC(self, rho, momentum, energy, B):
      rho[self.regionA.sliceBC] = rho[self.regionB.sliceOne]
      rho[self.regionB.sliceBC] = rho[self.regionA.sliceOne]
      
      
      if not par.staggered:
         momentum[0][self.regionA.sliceBC] = momentum[0][self.regionB.sliceOne]
         momentum[0][self.regionB.sliceBC] = momentum[0][self.regionA.sliceOne]

         if par.dim == 2:
            momentum[1][self.regionA.sliceBC] = momentum[1][self.regionB.sliceOne]
            momentum[1][self.regionB.sliceBC] = momentum[1][self.regionA.sliceOne]
            
      
      else: # Assuming that region A is 'R' or 'T'    HARDCODED
         if self.regionA.region=='R':
            if par.dim==1:
               momentum[0][-1] = momentum[0][self.regionB.sliceOne]
               momentum[0][self.regionB.sliceBC] = momentum[0][-2]
               
            else:
               momZL = momentum[0][ 1:-1, -2]
               momentum[0][ 1:-1, -2] = momentum[0][1:-1, 0 ]
               momentum[0][ 1:-1, 0] = momZL
               momentum[1][self.regionA.sliceBC] = momentum[1][self.regionB.sliceOne]
               momentum[1][self.regionB.sliceBC] = momentum[1][self.regionA.sliceOne]
               

         if self.regionA.region=='T':
            momentum[0][self.regionA.sliceBC] = momentum[0][self.regionB.sliceOne]
            momentum[0][self.regionB.sliceBC] = momentum[0][self.regionA.sliceOne]
            momYT = momentum[1][ -2,1:-1]
            momentum[1][-2,1:-1] = momentum[1][0,1:-1]
            momentum[1][0, 1:-1] = momYT   
            
            
      
      # Temporary placement, should be changed for each case (also Bz for B[0] and so on
      # For staggered, MHD:
      if self.regionA.region=='R':
         energy[1:-1,-1] = var.P[1:-1,1]/(par.gamma-1.) + 0.5*(var.Bz[1:-1,1]*var.Bz[1:-1,1] + var.By[1:-1,1]*var.By[1:-1,1] + var.Bx[1:-1,1]*var.Bx[1:-1,1])   \
         + 0.5*( momentum[0][1:-1,0]*momentum[0][1:-1,0] + 0.5*(momentum[1][1:-1,1]+momentum[1][:-2,1])*0.5*(momentum[1][1:-1,1]+momentum[1][:-2,1])  + momentum[2][1:-1,1]*momentum[2][1:-1,1] )/rho[1:-1,1]
         
         energy[1:-1, 0] = var.P[1:-1,-2]/(par.gamma-1.) + 0.5*(var.Bz[1:-1,-2]*var.Bz[1:-1,-2] + var.By[1:-1,-2]*var.By[1:-1,-2] + var.Bx[1:-1,-2]*var.Bx[1:-1,-2])   \
         + 0.5*( momentum[0][1:-1,-2]*momentum[0][1:-1,-2] + 0.5*(momentum[1][1:-1,-2]+momentum[1][:-2,-2])*0.5*(momentum[1][1:-1,-2]+momentum[1][:-2,-2])  + momentum[2][1:-1,-2]*momentum[2][1:-1,-2] )/rho[1:-1,-2]
         
      if self.regionA.region=='T':
         energy[-1,1:-1] = var.P[1,1:-1]/(par.gamma-1.) + 0.5*(var.Bz[1,1:-1]*var.Bz[1,1:-1] + var.By[1,1:-1]*var.By[1,1:-1] + var.Bx[1,1:-1]*var.Bx[1,1:-1])   \
         + 0.5*( 0.5*(momentum[0][1,1:-1]+momentum[0][1,:-2])*0.5*(momentum[0][1,1:-1]+momentum[0][1,:-2]) + momentum[1][0,1:-1]*momentum[1][0,1:-1]  + momentum[2][1,1:-1]*momentum[2][1,1:-1] )/rho[1,1:-1]
         
         energy[0,1:-1] = var.P[-2,1:-1]/(par.gamma-1.) + 0.5*(var.Bz[-2,1:-1]*var.Bz[-2,1:-1] + var.By[-2,1:-1]*var.By[-2,1:-1] + var.Bx[-2,1:-1]*var.Bx[-2,1:-1])   \
         + 0.5*( 0.5*(momentum[0][-2,1:-1]+momentum[0][-2,:-2])*0.5*(momentum[0][-2,1:-1]+momentum[0][-2,:-2]) + momentum[1][-2,1:-1]*momentum[1][-2,1:-1]  + momentum[2][-2,1:-1]*momentum[2][-2,1:-1] )/rho[-2,1:-1]         
   
      
      
      #energy[self.regionA.sliceBC] = energy[self.regionB.sliceOne]
      #energy[self.regionB.sliceBC] = energy[self.regionA.sliceOne]

      if par.MHD:
         momentum[2][self.regionA.sliceBC] = momentum[2][self.regionB.sliceOne]
         momentum[2][self.regionB.sliceBC] = momentum[2][self.regionA.sliceOne]

         B[0][self.regionA.sliceBC] = B[0][self.regionB.sliceOne]
         B[0][self.regionB.sliceBC] = B[0][self.regionA.sliceOne]
         
         B[1][self.regionA.sliceBC] = B[1][self.regionB.sliceOne]
         B[1][self.regionB.sliceBC] = B[1][self.regionA.sliceOne]

         B[2][self.regionA.sliceBC] = B[2][self.regionB.sliceOne]
         B[2][self.regionB.sliceBC] = B[2][self.regionA.sliceOne]


