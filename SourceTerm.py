#--------------------------------------------
#		Source.py
# This module computes the momentum and energy
# source term that should be added to the 
# conservative part of the equations. This can be
# further improved using Operator Splitting
#
#--------------------------------------------

import numpy as np
import Parameters as par
import Variables as var
import Settings as sets
import Grid
import ChangeOfVar
from scipy import linalg

print 'Loading SourceTerm..'

def computeGravSource():
   if par.GravityMode == 'Constant':
      momentumGravSourceZ = var.rho+var.lastrho*par.g
      momentumGravSourceY = var.rho*0.
      energyGravSource = var.momentumZ*par.g

   elif par.GravityMode == 'Radial':
      R2 = Grid.z*Grid.z + Grid.y*Grid.y + 0.005
      unitZ = Grid.z/np.sqrt(R2)
      unitY = Grid.y/np.sqrt(R2)
      g_Z = unitZ*par.g/R2 
      g_Y = unitY*par.g/R2
      momentumGravSourceZ = 0.5*(var.rho+var.lastrho)*g_Z
      momentumGravSourceY = 0.5*(var.rho+var.lastrho)*g_Y
      energyGravSource = var.momentumZ*g_Z + var.momentumY*g_Y
   return momentumGravSourceZ, momentumGravSourceY, energyGravSource




def computeMomentumDamping():

   momentumDampingSource = np.zeros(Grid.z.shape)
   energyDampingSource = np.zeros(Grid.z.shape)

   if par.DampingMode == 'Simple':
      DampingVel = par.DampingMultiplier*var.vZ   #Old method
    
      momentumDampingSource = -DampingVel*var.rho
      energyDampingSource = -0.5*DampingVel*DampingVel*var.rho

   
   elif par.DampingMode == 'FreeFallTime':

      H_p = var.P[:-1]/np.abs(var.P[1:]-var.P[:-1]) * Grid.dz  #pressure scale-height
      H_p = np.append(H_p, H_p[-1]) #Extent the array to avoid different sizes (this point is at the boundary)

      tau_ff = np.sqrt(2*H_p/np.abs(par.g))    #Time taken to fall 4Mm
      tau_damp = tau_ff

      n_iter = tau_damp[1:-1]/par.dt

      DampingPercent = 1./n_iter

      DampingPercent_scalar = np.max(DampingPercent)* par.DampingMultiplier

      DampingVel = DampingPercent_scalar*var.vZ


   elif par.DampingMode == 'MaximumVelocity':
      velMask = np.abs(var.vZ) > par.DampingMaxVel
      DampingVel = par.DampingMultiplier*var.vZ*velMask   #The damping is only applied for points with excess of velocity
 
      momentumDampingSource = -DampingVel*var.rho
      energyDampingSource = -0.5*DampingVel*DampingVel*var.rho


   else: 
      print '#### ERROR #### \t ' + par.DampingMode + ' is not defined as a DampingMode' 
   
   """ DELETE SOONER OR LATER
   H_p = var.P[:-1]/np.abs(var.P[1:]-var.P[:-1]) * Grid.dz  #pressure scale-height
   H_p = np.append(H_p, H_p[-1]) #Extent the array to avoid different sizes (this point is at the boundary)

   tau_ff = np.sqrt(2*H_p/np.abs(par.g))    #Time taken to fall 4Mm
   tau_damp = tau_ff
   
   n_iter = tau_damp[1:-1]/par.dt
   
   DampingPercent = 1./n_iter 
   
   DampingPercent_scalar = np.max(DampingPercent)* par.DampingMultiplier  
   
   if DampingPercent_scalar >1.:
      print 'Too much Damping!!'
   
   DampingVel = DampingPercent_scalar*var.v
   
   #print np.argmax(DampingPercent)
   #DampingVel = par.DampingMultiplier*var.v   #Old method
   
   #momentumDampingSource = -DampingVel*var.rho
   #energyDampingSource = -0.5*DampingVel*DampingVel*var.rho
   """
   return momentumDampingSource, energyDampingSource

def computeTemperatureDiffusion():
   #ct = 9e-12 * 1e5 # Value at E.Priest "Solar Magnetohydrodynamics" * mks to cgs factor
   if par.SpitzerDiffusion:
      var.kappa = par.ct*var.T**(5./2.)
   
   Dkappa = (var.kappa[2:]-var.kappa[:-2])/4.
   EnergyDiff = (var.kappa[1:-1]+Dkappa)*var.T[2:] - 2.*var.kappa[1:-1]*var.T[1:-1] + (var.kappa[1:-1]-Dkappa)*var.T[:-2]

   return par.DiffusionPercent*EnergyDiff/(Grid.dz*Grid.dz)



def computeImplicitConduction():
   if par.SpitzerDiffusion:
      var.kappa = par.ct*var.T**(5./2.)
      
   """   
   #Hard-coded boundaries
   var.rho[0] = 2.*var.rho[1]-var.rho[2]
   boundaryRho = 0.5*(var.rho[0]+var.rho[1])
   var.momentum[0] = -var.momentum[1]
   
   
   var.rho[-1] = var.rho[-2]
   var.momentum[-1] = var.momentum[-2]
   var.v = var.momentum/var.rho
   """ 
      
   e = par.cv*var.T   #Internal energy
   
   Dkappa = (var.kappa[2:]-var.kappa[:-2])/4.
   
   dz2 = Grid.dz*Grid.dz
   A = par.DiffusionPercent * (var.kappa[1:-1]+Dkappa)/dz2
   B = -2.*par.DiffusionPercent * var.kappa[1:-1]/dz2
   C = par.DiffusionPercent * (var.kappa[1:-1]-Dkappa)/dz2
   
   var.diag[1:-1] = par.dt*B/(var.rho[1:-1]*par.cv) - 1.  
   var.lower[:-2] = par.dt*C/(var.rho[:-2]*par.cv)   
   var.upper[2:] = par.dt*A/(var.rho[2:]*par.cv)
   
   var.rhs[1:-1] = -e[1:-1]*var.rho[1:-1]
   
   if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC()
   if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC()
   
   """
   # Hard-coded fixed temperature boundaries
   #Lower boundary
   diag = np.insert(diag, 0, 1.)   
   upper = np.insert(upper, 0, 1.)
   meanrho = 0.5*(var.rho[0] + var.rho[1])
   rhs[0] = 2.*meanrho*par.cv*1e4
   
   #Upper boundary
   diag = np.append(diag, 1.)
   lower = np.append(lower, 1.)
   meanrho = 0.5*(var.rho[-2] + var.rho[-1])
   rhs[-1] = 2.*meanrho*par.cv*1e6
   
   #Fill of the lower and upper diagonals https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html#scipy.linalg.solve_banded
   #upper = np.insert(upper, 0 , 0.)
   #lower = np.append(lower, 0.)
   """
   #print e[:6]*var.rho[:6]
   sol_e = linalg.solve_banded( (1,1), [var.upper, var.diag, var.lower], var.rhs )
   #print sol_e[:6]
   #print ''
   #print var.T-sol_e/(par.cv*var.rho)
   if par.dim == 1:
      var.energy = sol_e + 0.5*var.momentumZ*var.momentumZ/var.rho
   elif par.dim == 2:
      var.energy = sol_e + 0.5*( var.momentumZ*var.momentumZ + var.momentumY*var.momentumY )/var.rho
   
   var.T = sol_e/(var.rho*par.cv)


def computeRadiativeLosses():
   logT = np.log10(var.T)

   """
   Lamda = np.zeros(Grid.z.shape)  #This can be allocated at settings/variables
   for i in range(len(logT)):   #Not sure if this will be faster than masked arrays...
     if logT[i] <= 4.97:
       Lamda[i] = 1.09e-31*var.T[i]*var.T[i]
     elif logT[i] < 5.67:
       Lamda[i] = 8.87e-17/var.T[i]
     elif logT[i] < 6.18:
       Lamda[i] = 1.90e-22
     elif logT[i] < 6.55:
       Lamda[i] = 3.53e-13*var.T[i]**(-3./2.)
     elif logT[i] < 6.90:
       Lamda[i] = 3.56e-25*var.T[i]**(1./3.)
     else:
       Lamda[i] = 5.49e-16/var.T[i]
   logLamda = np.log10(Lamda)
   """
   
   logLamda = np.interp(logT, var.logT_table, var.logLamda_table)
   lowTmask = logT<4.4771212547196626  #log10(3e4)
   logLamda[lowTmask] = var.logLamda_table[9] + 35*(logT[lowTmask]- 4.4771212547196626)

   numericalDensity = var.rho/(par.mu*par.molarMass)
   var.RadiativeLoss = par.RadiationPercent * numericalDensity * numericalDensity * 10**logLamda
   #var.RadiativeLoss = var.RadiativeLoss*(var.T>3e4)
   return var.RadiativeLoss


def ComputeSource():

   if par.IsThereGravity:
      momentumGZ, momentumGY, energyG = computeGravSource()
      var.momentumZ += par.dt*momentumGZ
      var.momentumY += par.dt*momentumGY
      var.energy += par.dt*energyG
      ChangeOfVar.ConvertToPrim()


  
   if par.RadiativeLoss:
      RadLoss = computeRadiativeLosses()
      var.energy[1:-1] -= par.dt*RadLoss[1:-1]
      ChangeOfVar.ConvertToPrim()

   if par.ThermalDiffusion:
      if par.ImplicitConduction:
         computeImplicitConduction()
      else:
         ThermalDiff = computeTemperatureDiffusion()
         var.energy[1:-1] += par.dt*ThermalDiff
      ChangeOfVar.ConvertToPrim()
       
   if par.MomentumDamping:
      momentumDamping, energyDamping = computeMomentumDamping()
      var.momentum += momentumDamping
      var.energy += energyDamping
      ChangeOfVar.ConvertToPrim()
     
     
     
     
     
     
     
