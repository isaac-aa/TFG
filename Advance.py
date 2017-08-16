#--------------------------------------------
#		Advance.py
# This module defines the different numerical
# schemes that can be used in this code. The
# scheme should be set at Settings.py
#
#--------------------------------------------



import Parameters as par
import Settings as sets 
import Grid
import Variables as var
import numpy as np
import Flux
import SourceTerm
import ChangeOfVar



class AdvanceScheme():
   """
   Base class for the numerical scheme for the advance of
   the conservative part of the equations, ie. without source terms,
   gravity and thermal conduction
   """
   
   name = 'General Scheme'
   
   def __str__(self):
      return 'Numerical scheme: ' + self.name
   
   
   def setup(self):
      """
      Function to pre-allocate all the arrays, if needed
      """
      print 'Setup \t ' + str(self)

   def compute(self):
      """
      Main implementation of the scheme.
      Need to be overriden for each scheme
      """
      print '##WARNING## \n \n Select a BC'
      



class LaxFriedrichs1D(AdvanceScheme):
   name = 'LaxFriedrichs1D'
        

   def compute(self):
      self.lamdaZ = par.dt/Grid.dz
      
      massFlux, momentumFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      var.rho[1:-1] = 0.5*(var.rho[2:]+var.rho[:-2]) - 0.5*self.lamdaZ*(massFlux[2:] - massFlux[:-2]) 
      var.momentumZ[1:-1] = 0.5*(var.momentumZ[2:]+var.momentumZ[:-2]) - 0.5*self.lamdaZ*(momentumFlux[2:] - momentumFlux[:-2]) 
      var.energy[1:-1] = 0.5*(var.energy[2:]+var.energy[:-2]) - 0.5*self.lamdaZ*(energyFlux[2:] - energyFlux[:-2]) 



class LaxFriedrichs2D(AdvanceScheme):
   name = 'LaxFriedrichs2D'
   

   def compute(self):
      self.lamdaZ = par.dt/Grid.dz
      self.lamdaY = par.dt/Grid.dy
      
      massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY = Flux.ComputeFlux2D(var.rho, var.momentumZ, var.momentumY, var.energy)

      var.rho[1:-1, 1:-1] = 0.25*(var.rho[2:, 1:-1] + var.rho[:-2, 1:-1] + var.rho[1:-1, 2:] + var.rho[1:-1, :-2]) - 0.5*self.lamdaZ*(massFluxZ[1:-1, 2:] - massFluxZ[1:-1, :-2]) - 0.5*self.lamdaY*(massFluxY[2:, 1:-1] - massFluxY[:-2, 1:-1])
      var.momentumZ[1:-1, 1:-1] = 0.25*(var.momentumZ[2:, 1:-1] + var.momentumZ[:-2, 1:-1] + var.momentumZ[1:-1, 2:] + var.momentumZ[1:-1, :-2]) - 0.5*self.lamdaZ*(momentumZFluxZ[1:-1, 2:] - momentumZFluxZ[1:-1, :-2]) - 0.5*self.lamdaY*(momentumZFluxY[2:, 1:-1] - momentumZFluxY[:-2, 1:-1])
      var.momentumY[1:-1, 1:-1] = 0.25*(var.momentumY[2:, 1:-1] + var.momentumY[:-2, 1:-1] + var.momentumY[1:-1, 2:] + var.momentumY[1:-1, :-2]) - 0.5*self.lamdaY*(momentumYFluxZ[1:-1, 2:] - momentumYFluxZ[1:-1, :-2]) - 0.5*self.lamdaY*(momentumYFluxY[2:, 1:-1] - momentumYFluxY[:-2, 1:-1]) 
      var.energy[1:-1, 1:-1] = 0.25*(var.energy[2:, 1:-1] + var.energy[:-2, 1:-1] + var.energy[1:-1, 2:] + var.energy[1:-1, :-2]) - 0.5*self.lamdaZ*(energyFluxZ[1:-1, 2:] - energyFluxZ[1:-1, :-2]) - 0.5*self.lamdaY*(energyFluxY[2:, 1:-1] - energyFluxY[:-2, 1:-1])



class FirstGen(AdvanceScheme):
   name = 'Lax-Wendroff Ritchmyer 1D'
   
   def setup(self):
      print 'Setup \t ' + str(self)
      
      self.rhoHalf = np.ones_like(Grid.z)
      self.momentumHalf = np.ones_like(Grid.z)
      self.energyHalf = np.ones_like(Grid.z)
      
      self.vHalf = np.ones_like(Grid.z)
      self.PHalf = np.ones_like(Grid.z)


   def compute(self):
      
      self.lamdaZ = par.dt/Grid.dz

      massFlux, momentumFlux, energyFlux = Flux.ComputeFluxOld(var.rho, var.momentumZ, var.energy, var.vZ, var.P)

      self.rhoHalf[:-1] = 0.5*(var.rho[1:]+var.rho[:-1]) - 0.5*self.lamdaZ*(massFlux[1:] - massFlux[:-1]) 
      self.momentumHalf[:-1] = 0.5*(var.momentumZ[1:]+var.momentumZ[:-1]) - 0.5*self.lamdaZ*(momentumFlux[1:] - momentumFlux[:-1]) 
      self.energyHalf[:-1] = 0.5*(var.energy[1:]+var.energy[:-1]) - 0.5*self.lamdaZ*(energyFlux[1:] - energyFlux[:-1]) 


      vHalf = self.momentumHalf/self.rhoHalf
      eHalf = self.energyHalf/self.rhoHalf - 0.5*vHalf*vHalf              #rho*e = E - 0.5*rho*v*v
      PHalf = self.rhoHalf*(par.gamma-1.)*eHalf
      
      massFluxHalf, momentumFluxHalf, energyFluxHalf = Flux.ComputeFluxOld(self.rhoHalf, self.momentumHalf, self.energyHalf, vHalf, PHalf)


      var.rho[1:-1] = var.rho[1:-1] - self.lamdaZ*( massFluxHalf[1:-1] - massFluxHalf[:-2] )
      var.momentumZ[1:-1] = var.momentumZ[1:-1] - self.lamdaZ*( momentumFluxHalf[1:-1] - momentumFluxHalf[:-2] ) 
      var.energy[1:-1] = var.energy[1:-1] - self.lamdaZ*( energyFluxHalf[1:-1] - energyFluxHalf[:-2] )



class RK4(AdvanceScheme):
   name = 'Runge-Kutta 4th order'
   
   def setup(self):
      self.rho_k1 = np.zeros_like(Grid.z)
      self.rho_k2 = np.zeros_like(Grid.z)
      self.rho_k3 = np.zeros_like(Grid.z)
      self.rho_k4 = np.zeros_like(Grid.z)
      
      self.momZ_k1 = np.zeros_like(Grid.z)
      self.momZ_k2 = np.zeros_like(Grid.z)
      self.momZ_k3 = np.zeros_like(Grid.z)
      self.momZ_k4 = np.zeros_like(Grid.z)
      
      self.ene_k1 = np.zeros_like(Grid.z)
      self.ene_k2 = np.zeros_like(Grid.z)
      self.ene_k3 = np.zeros_like(Grid.z)
      self.ene_k4 = np.zeros_like(Grid.z)
      
      if par.dim == 2:
         self.momY_k1 = np.zeros_like(Grid.z)
         self.momY_k2 = np.zeros_like(Grid.z)
         self.momY_k3 = np.zeros_like(Grid.z)
         self.momY_k4 = np.zeros_like(Grid.z)
      
      AdvanceScheme.setup(self)
   
   def compute(self):
      if par.dim==1:
         self.compute1D()
      elif par.dim==2:
         self.compute2D()
   

   def compute1D(self):
      h = par.dt


      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k1[1:-1] = -(massFlux[2:]-massFlux[:-2])/(2.*Grid.dz) 
      self.momZ_k1[1:-1] = -(momentumZFlux[2:]-momentumZFlux[:-2])/(2.*Grid.dz) + momentumGZ[1:-1]
      self.ene_k1[1:-1] = -(energyFlux[2:] - energyFlux[:-2])/(2.*Grid.dz) - SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      
      var.rho = var.lastrho+self.rho_k1*h*0.5
      var.momentumZ = var.lastmomentumZ+self.momZ_k1*h*0.5
      var.energy = var.lastenergy+self.ene_k1*h*0.5
      ChangeOfVar.ConvertToPrim()
      
      if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
      if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)

      
      
      
      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)
      
      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k2[1:-1] = -(massFlux[2:]-massFlux[:-2])/(2.*Grid.dz) 
      self.momZ_k2[1:-1] = -(momentumZFlux[2:]-momentumZFlux[:-2])/(2.*Grid.dz) + momentumGZ[1:-1]
      self.ene_k2[1:-1] = -(energyFlux[2:] - energyFlux[:-2])/(2.*Grid.dz) - SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      
      var.rho = var.lastrho+self.rho_k2*h*0.5
      var.momentumZ = var.lastmomentumZ+self.momZ_k2*h*0.5
      var.energy = var.lastenergy+self.ene_k2*h*0.5
      ChangeOfVar.ConvertToPrim()
      
      if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
      if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)



      
      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k3[1:-1] = -(massFlux[2:]-massFlux[:-2])/(2.*Grid.dz) 
      self.momZ_k3[1:-1] = -(momentumZFlux[2:]-momentumZFlux[:-2])/(2.*Grid.dz) + momentumGZ[1:-1]
      self.ene_k3[1:-1] = -(energyFlux[2:] - energyFlux[:-2])/(2.*Grid.dz) - SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      
      var.rho = var.lastrho+self.rho_k3*h*0.5
      var.momentumZ = var.lastmomentumZ+self.momZ_k3*h*0.5
      var.energy = var.lastenergy+self.ene_k3*h*0.5
      ChangeOfVar.ConvertToPrim()
      
      if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
      if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)

      
      
      
      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k4[1:-1] = -(massFlux[2:]-massFlux[:-2])/(2.*Grid.dz) 
      self.momZ_k4[1:-1] = -(momentumZFlux[2:]-momentumZFlux[:-2])/(2.*Grid.dz) + momentumGZ[1:-1]
      self.ene_k4[1:-1] = -(energyFlux[2:] - energyFlux[:-2])/(2.*Grid.dz) - SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]


      var.rho = var.lastrho +              h/6. * (self.rho_k1 + 2.*self.rho_k2 + 2.*self.rho_k3 + self.rho_k4)
      var.momentumZ = var.lastmomentumZ +  h/6. * (self.momZ_k1 + 2.*self.momZ_k2 + 2.*self.momZ_k3 + self.momZ_k4)
      var.energy = var.lastenergy +        h/6. * (self.ene_k1 + 2.*self.ene_k2 + 2.*self.ene_k3 + self.ene_k4)



   def compute2D(self):
      h = par.dt
      
      massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY = Flux.ComputeFlux2D(var.rho, var.momentumZ, var.momentumY, var.energy)


      self.rho_k1[1:-1,1:-1] = -(massFluxZ[1:-1, 2:]-massFluxZ[1:-1, :-2])/(2.*Grid.dz) - (massFluxY[2:, 1:-1]-massFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.momZ_k1[1:-1,1:-1] = -(momentumZFluxZ[1:-1, 2:]-momentumZFluxZ[1:-1, :-2])/(2.*Grid.dz) - (momentumZFluxY[2:,1:-1]-momentumZFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.momY_k1[1:-1,1:-1] = -(momentumYFluxZ[1:-1, 2:]-momentumYFluxZ[1:-1, :-2])/(2.*Grid.dz) - (momentumYFluxY[2:,1:-1]-momentumYFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.ene_k1[1:-1,1:-1] = -(energyFluxZ[1:-1, 2:] - energyFluxZ[1:-1, :-2])/(2.*Grid.dz) - (energyFluxY[2:,1:-1] - energyFluxY[:-2,1:-1])/(2.*Grid.dy)

      massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY = Flux.ComputeFlux2D(var.rho+self.rho_k1*h*0.5, var.momentumZ+self.momZ_k1*h*0.5, var.momentumY+self.momY_k1*h*0.5, var.energy+self.ene_k1*h*0.5)

      self.rho_k2[1:-1,1:-1] = -(massFluxZ[1:-1, 2:]-massFluxZ[1:-1, :-2])/(2.*Grid.dz) - (massFluxY[2:, 1:-1]-massFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.momZ_k2[1:-1,1:-1] = -(momentumZFluxZ[1:-1, 2:]-momentumZFluxZ[1:-1, :-2])/(2.*Grid.dz) - (momentumZFluxY[2:,1:-1]-momentumZFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.momY_k2[1:-1,1:-1] = -(momentumYFluxZ[1:-1, 2:]-momentumYFluxZ[1:-1, :-2])/(2.*Grid.dz) - (momentumYFluxY[2:,1:-1]-momentumYFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.ene_k2[1:-1,1:-1] = -(energyFluxZ[1:-1, 2:] - energyFluxZ[1:-1, :-2])/(2.*Grid.dz) - (energyFluxY[2:,1:-1] - energyFluxY[:-2,1:-1])/(2.*Grid.dy)


      massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY = Flux.ComputeFlux2D(var.rho+self.rho_k2*h*0.5, var.momentumZ+self.momZ_k2*h*0.5, var.momentumY+self.momY_k2*h*0.5, var.energy+self.ene_k2*h*0.5)

      self.rho_k3[1:-1,1:-1] = -(massFluxZ[1:-1, 2:]-massFluxZ[1:-1, :-2])/(2.*Grid.dz) - (massFluxY[2:, 1:-1]-massFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.momZ_k3[1:-1,1:-1] = -(momentumZFluxZ[1:-1, 2:]-momentumZFluxZ[1:-1, :-2])/(2.*Grid.dz) - (momentumZFluxY[2:,1:-1]-momentumZFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.momY_k3[1:-1,1:-1] = -(momentumYFluxZ[1:-1, 2:]-momentumYFluxZ[1:-1, :-2])/(2.*Grid.dz) - (momentumYFluxY[2:,1:-1]-momentumYFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.ene_k3[1:-1,1:-1] = -(energyFluxZ[1:-1, 2:] - energyFluxZ[1:-1, :-2])/(2.*Grid.dz) - (energyFluxY[2:,1:-1] - energyFluxY[:-2,1:-1])/(2.*Grid.dy)


      massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY = Flux.ComputeFlux2D(var.rho+self.rho_k3*h, var.momentumZ+self.momZ_k3*h, var.momentumY+self.momY_k3*h, var.energy+self.ene_k3*h)

      self.rho_k4[1:-1,1:-1] = -(massFluxZ[1:-1, 2:]-massFluxZ[1:-1, :-2])/(2.*Grid.dz) - (massFluxY[2:, 1:-1]-massFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.momZ_k4[1:-1,1:-1] = -(momentumZFluxZ[1:-1, 2:]-momentumZFluxZ[1:-1, :-2])/(2.*Grid.dz) - (momentumZFluxY[2:,1:-1]-momentumZFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.momY_k4[1:-1,1:-1] = -(momentumYFluxZ[1:-1, 2:]-momentumYFluxZ[1:-1, :-2])/(2.*Grid.dz) - (momentumYFluxY[2:,1:-1]-momentumYFluxY[:-2,1:-1])/(2.*Grid.dy)
      self.ene_k4[1:-1,1:-1] = -(energyFluxZ[1:-1, 2:] - energyFluxZ[1:-1, :-2])/(2.*Grid.dz) - (energyFluxY[2:,1:-1] - energyFluxY[:-2,1:-1])/(2.*Grid.dy)


      var.rho = var.rho +              h/6. * (self.rho_k1 + 2*self.rho_k2 + 2*self.rho_k3 + self.rho_k4)
      var.momentumZ = var.momentumZ +  h/6. * (self.momZ_k1 + 2*self.momZ_k2 + 2*self.momZ_k3 + self.momZ_k4)
      var.momentumY = var.momentumY +  h/6. * (self.momY_k1 + 2*self.momY_k2 + 2*self.momY_k3 + self.momY_k4)
      var.energy = var.energy +        h/6. * (self.ene_k1 + 2*self.ene_k2 + 2*self.ene_k3 + self.ene_k4)








class RK3(AdvanceScheme):
   name = 'Runge-Kutta 3th order'
   
   def setup(self):
      self.rho_k1 = np.zeros_like(Grid.z)
      self.rho_k2 = np.zeros_like(Grid.z)
      self.rho_k3 = np.zeros_like(Grid.z)
      
      self.momZ_k1 = np.zeros_like(Grid.z)
      self.momZ_k2 = np.zeros_like(Grid.z)
      self.momZ_k3 = np.zeros_like(Grid.z)
      
      self.ene_k1 = np.zeros_like(Grid.z)
      self.ene_k2 = np.zeros_like(Grid.z)
      self.ene_k3 = np.zeros_like(Grid.z)
      
      if par.dim == 2:
         self.momY_k1 = np.zeros_like(Grid.z)
         self.momY_k2 = np.zeros_like(Grid.z)
         self.momY_k3 = np.zeros_like(Grid.z)
      
      AdvanceScheme.setup(self)
   
   def compute(self):
      if par.dim==1:
         self.compute1D()
   

   def compute1D(self):
      h = par.dt

      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k1[1:-1] = -(massFlux[2:]-massFlux[:-2])/(2.*Grid.dz) 
      self.momZ_k1[1:-1] = -(momentumZFlux[2:]-momentumZFlux[:-2])/(2.*Grid.dz) #+ momentumGZ[1:-1]
      self.ene_k1[1:-1] = -(energyFlux[2:] - energyFlux[:-2])/(2.*Grid.dz) #- SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      
      var.rho = var.lastrho+self.rho_k1*h*0.5
      var.momentumZ = var.lastmomentumZ+self.momZ_k1*h*0.5
      var.energy = var.lastenergy+self.ene_k1*h*0.5
      ChangeOfVar.ConvertToPrim()
      
      if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
      if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)

      
      
      
      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)
      
      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k2[1:-1] = -(massFlux[2:]-massFlux[:-2])/(2.*Grid.dz) 
      self.momZ_k2[1:-1] = -(momentumZFlux[2:]-momentumZFlux[:-2])/(2.*Grid.dz) #+ momentumGZ[1:-1]
      self.ene_k2[1:-1] = -(energyFlux[2:] - energyFlux[:-2])/(2.*Grid.dz) #- SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      
      var.rho = var.lastrho - self.rho_k1*h + 2.*self.rho_k2*h
      var.momentumZ = var.lastmomentumZ - self.momZ_k1*h + 2.*self.momZ_k2*h
      var.energy = var.lastenergy - self.ene_k1*h + 2.*self.ene_k2*h
      ChangeOfVar.ConvertToPrim()
      
      if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
      if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)



      
      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k3[1:-1] = -(massFlux[2:]-massFlux[:-2])/(2.*Grid.dz) 
      self.momZ_k3[1:-1] = -(momentumZFlux[2:]-momentumZFlux[:-2])/(2.*Grid.dz) #+ momentumGZ[1:-1]
      self.ene_k3[1:-1] = -(energyFlux[2:] - energyFlux[:-2])/(2.*Grid.dz) #- SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      

      var.rho = var.lastrho +              h/6. * (self.rho_k1 + 4.*self.rho_k2 + self.rho_k3)
      var.momentumZ = var.lastmomentumZ +  h/6. * (self.momZ_k1 + 4.*self.momZ_k2 + self.momZ_k3)
      var.energy = var.lastenergy +        h/6. * (self.ene_k1 + 4.*self.ene_k2 + self.ene_k3)





class RK3Staggered(AdvanceScheme):
   name = 'Runge-Kutta 3th order on staggered grid'
   
   def setup(self):
      self.rho_k1 = np.zeros_like(Grid.z)
      self.rho_k2 = np.zeros_like(Grid.z)
      self.rho_k3 = np.zeros_like(Grid.z)
      
      self.momZ_k1 = np.zeros_like(Grid.z)
      self.momZ_k2 = np.zeros_like(Grid.z)
      self.momZ_k3 = np.zeros_like(Grid.z)
      
      self.ene_k1 = np.zeros_like(Grid.z)
      self.ene_k2 = np.zeros_like(Grid.z)
      self.ene_k3 = np.zeros_like(Grid.z)
      
      if par.dim == 2:
         self.momY_k1 = np.zeros_like(Grid.z)
         self.momY_k2 = np.zeros_like(Grid.z)
         self.momY_k3 = np.zeros_like(Grid.z)
      
      AdvanceScheme.setup(self)
   
   def compute(self):
      if par.dim==1:
         self.compute1D()
   

   def compute1D(self):
      h = par.dt

      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k1[1:-1] = -(var.momentumZ[1:-1]-var.momentumZ[:-2])/Grid.dz
      
      term1 = (0.5*(var.momentumZ[:-2] + var.momentumZ[1:-1]))**2/var.rho[1:-1]
      self.momZ_k1[1:-2] = -(term1[1:]-term1[:-1])/Grid.dz - (var.P[2:-1]-var.P [1:-2])/Grid.dz
      
      Et_plus_P = (0.5*(var.energy+var.P)[:-1] + (var.energy+var.P)[1:])/(0.5*(var.rho[1:] + var.rho[:-1])  ) 
      term = Et_plus_P*var.momentumZ[:-1]
      self.ene_k1[1:-1] = -(term[1:]-term[:-1])/Grid.dz #- SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      
      var.rho = var.lastrho+self.rho_k1*h*0.5
      var.momentumZ = var.lastmomentumZ+self.momZ_k1*h*0.5
      var.energy = var.lastenergy+self.ene_k1*h*0.5
      ChangeOfVar.ConvertToPrim()
      
      if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
      if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)

      
      
      
      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)
      
      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k2[1:-1] = -(var.momentumZ[1:-1]-var.momentumZ[:-2])/Grid.dz
      term1 = (0.5*(var.momentumZ[:-2] + var.momentumZ[1:-1]))**2/var.rho[1:-1]
      self.momZ_k2[1:-2] = -(term1[1:]-term1[:-1])/Grid.dz - (var.P[2:-1]-var.P [1:-2])/Grid.dz
      Et_plus_P = (0.5*(var.energy+var.P)[:-1] + (var.energy+var.P)[1:])/(0.5*(var.rho[1:] + var.rho[:-1])  ) 
      term = Et_plus_P*var.momentumZ[0:-1]
      self.ene_k2[1:-1] = -(term[1:]-term[:-1])/Grid.dz #- SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      
      var.rho = var.lastrho - self.rho_k1*h + 2.*self.rho_k2*h
      var.momentumZ = var.lastmomentumZ - self.momZ_k1*h + 2.*self.momZ_k2*h
      var.energy = var.lastenergy - self.ene_k1*h + 2.*self.ene_k2*h
      ChangeOfVar.ConvertToPrim()
      
      if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
      if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)



      
      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      self.rho_k3[1:-1] = -(var.momentumZ[1:-1]-var.momentumZ[:-2])/Grid.dz
      term1 = (0.5*(var.momentumZ[:-2] + var.momentumZ[1:-1]))**2/var.rho[1:-1]
      self.momZ_k3[1:-2] = -(term1[1:]-term1[:-1])/Grid.dz - (var.P[2:-1]-var.P [1:-2])/Grid.dz
      Et_plus_P = (0.5*(var.energy+var.P)[:-1] + (var.energy+var.P)[1:])/(0.5*(var.rho[1:] + var.rho[:-1])  ) 
      term = Et_plus_P*var.momentumZ[0:-1]
      self.ene_k3[1:-1] = -(term[1:]-term[:-1])/Grid.dz #- SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1]
      

      var.rho = var.lastrho +              h/6. * (self.rho_k1 + 4.*self.rho_k2 + self.rho_k3)
      var.momentumZ = var.lastmomentumZ +  h/6. * (self.momZ_k1 + 4.*self.momZ_k2 + self.momZ_k3)
      var.energy = var.lastenergy +        h/6. * (self.ene_k1 + 4.*self.ene_k2 + self.ene_k3)








class ForwardEuler(AdvanceScheme):
   name = 'Forward Euler'
   
   def compute(self):
      if par.dim==1:
         self.compute1D()
   

   def compute1D(self):
      massFlux, momentumZFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentumZ, var.energy)

      momentumGZ, momentumGY, energyG = SourceTerm.computeGravSource()
      var.rho[1:-1]        = var.lastrho[1:-1]         + par.dt*( -(massFlux[2:]-massFlux[:-2])/(2.*Grid.dz) ) 
      var.momentumZ[1:-1]  = var.lastmomentumZ[1:-1]   + par.dt*( -(momentumZFlux[2:]-momentumZFlux[:-2])/(2.*Grid.dz) + momentumGZ[1:-1] )
      var.energy[1:-1]     = var.lastenergy[1:-1]      + par.dt*( -(energyFlux[2:] - energyFlux[:-2])/(2.*Grid.dz) - SourceTerm.computeRadiativeLosses()[1:-1] + energyG[1:-1])
      




   
