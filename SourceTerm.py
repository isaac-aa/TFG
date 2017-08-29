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
import scipy.sparse.linalg

print 'Loading SourceTerm..'


def BuildSystemMatrix(center, T, B, L, R):
   N = Grid.z.shape[0]*Grid.z.shape[1]
   y_side = Grid.z.shape[0]
   x_side = Grid.z.shape[1]
   

   nnz = np.zeros(5*N)  #Four neighbours for each cell plus the cell itself
   col = np.zeros(5*N, dtype=int)
   offset = np.zeros(N+1, dtype=int)
   i=1

   for xi in range(x_side):
      for yj in range(y_side):
         #index = (yj*y_side + xi)*5
         cell = yj*x_side + xi
         index = cell*5
         #print cell, xi, yj
         nnz[index] = center[yj,xi] #Center cell
         col[index] = cell


         nnz[index + 1] = L[yj,xi] #Left
         col[index + 1] = cell-1  #Left neighbor


         nnz[index + 2] = R[yj,xi] #Right
         col[index + 2] = cell+1  #Right


         nnz[index + 3] = T[yj,xi] #Top
         col[index + 3] = cell+x_side #Top


         nnz[index + 4] = B[yj,xi] #Bottom
         col[index + 4] = cell-x_side


         offset[i] = offset[i-1]+5  #Improve performance-wise later on
         i+=1

   """ Periodic BC: sounds good, doesnt work (properly)
   for xi in range(len(x)):
     index = (y_side*0 + xi)*5  #Lower boundary
     nnz[index + 3] = TB       #Posible error con top/bottom indices
     col[index + 3] = N-(x_side-xi)

     index = (y_side*(y_side-1) + xi)*5 #Upper boundary
     nnz[index + 4] = TB
     col[index + 4] = xi


   for yj in range(len(y)):
     index = (y_side*yj + 0)*5 #Left boundary
     cell = y_side*yj
     index = cell*5

     nnz[index + 1] = LR
     col[index + 1] = cell + (x_side-1)

     index = (y_side*yj + (x_side-1))*5 #Right boundary
     cell = y_side*yj + (x_side-1)
     nnz[index + 2] = LR
     col[index + 2] = cell - (x_side-1)
   """

   # Zero-derivative potential boundary condition
   rhs, nnzC, nnzBC = sets.BoundaryConditionB.energyBC.computeImplicitBC()
   for xi in range(x_side):
      index = (x_side*0 + xi)*5     # B
      nnz[index] = nnzC
      nnz[index + 3] = nnzBC
      var.rhs[index/5] = rhs[xi]
      nnz[index + 4] = 0
      col[index + 4] = 0

   rhs, nnzC, nnzBC = sets.BoundaryConditionT.energyBC.computeImplicitBC()
   for xi in range(x_side):
      index = (x_side*(y_side-1) + xi)*5   # T
      nnz[index] = nnzC
      nnz[index + 4] = nnzBC
      var.rhs[index/5] = rhs[xi]
      nnz[index + 3] = 0
      col[index + 3] = 0

   rhs, nnzC, nnzBC = sets.BoundaryConditionL.energyBC.computeImplicitBC()
   for yj in range(1,y_side-1):
      index = (x_side*yj + 0)*5            # L
      nnz[index] = nnzC
      nnz[index+2] = nnzBC
      var.rhs[index/5] = rhs[yj]
      nnz[index + 1] = 0.
      col[index + 1] = 0
      
      
   rhs, nnzC, nnzBC = sets.BoundaryConditionR.energyBC.computeImplicitBC()
   for yj in range(1,y_side-1):
      index = (x_side*yj + (x_side-1))*5    # R
      nnz[index] = nnzC
      nnz[index + 1] = nnzBC
      var.rhs[index/5] = rhs[yj]
      nnz[index + 2] = 0.
      col[index + 2] = 0
      
   
   #print nnz.shape
   #print nnz[2000:2020]
   #print col[2000:2020]

   return nnz, col, offset




def computeGravSource():
   if par.GravityMode == 'Constant':
      if par.staggered:
         momentumGravSourceZ = 0.5*(var.rho[1:-1,:-1] + var.rho[1:-1,1:])*par.g
         momentumGravSourceY = 0.
         energyGravSource = 0.5*(var.momentumZ[1:-1, :-2] + var.momentumZ[1:-1, 1:-1])*par.g
      else:
         momentumGravSourceZ = 0.5*(var.rho+var.lastrho)*par.g
         momentumGravSourceY = var.rho*0.
         energyGravSource = 0.5*(var.momentumZ + var.lastmomentumZ)*par.g
      
   elif par.GravityMode == 'Radial':
      R2 = Grid.z*Grid.z + Grid.y*Grid.y + 0.005
      unitZ = Grid.z/np.sqrt(R2)
      unitY = Grid.y/np.sqrt(R2)
      g_Z = unitZ*par.g/R2 
      g_Y = unitY*par.g/R2
      momentumGravSourceZ = 0.5*(var.rho+var.lastrho)*g_Z
      momentumGravSourceY = 0.5*(var.rho+var.lastrho)*g_Y
      energyGravSource = var.momentumZ*g_Z + var.momentumY*g_Y

   elif par.GravityMode == 'Poisson':
      """
      Esto se puede mejorar mucho usando:
      OOP, Cython, o al menos otra funcion
      """

      phi = np.zeros(Grid.z.shape) #NO REALOCAR ESTO!!!!
      fz = np.zeros(Grid.z.shape)
      fy = np.zeros(Grid.z.shape) 

      # Phi potential is computed:
      for i in np.arange(par.Nz-1, 0, -1):
         for j in np.arange(par.Ny-1, 0, -1):

            phi_iplus1 = 0.
            phi_jplus1 = 0.

            if i>3 and i<par.Nz-3 and j>3 and j<par.Ny-3:
               phi_iplus1 = ( var.rho[i-1,j] - (phi[i-1,j-1] - 2.*phi[i-1,j] + phi[i-1, j+1] )/(Grid.dy*Grid.dy) )*Grid.dz*Grid.dz - phi[i-2, j] + 2.*phi[i-1, j]
               phi_jplus1 = ( var.rho[i,j-1] - (phi[i-2,j-1] - 2.*phi[i-1,j-1] + phi[i, j-1] )/(Grid.dz*Grid.dz) )*Grid.dy*Grid.dy - phi[i, j-2] + 2.*phi[i, j-1]

            else:
               phi_iplus1 = 0. #( rho[i-1,j] )*dx*dx
               phi_jplus1 = 0. #( rho[i,j-1] )*dy*dy


            phi[i, j] = 0.5*(phi_iplus1+phi_jplus1)

      fz[1:-1, 1:-1] = (phi[2:, 1:-1] - phi[:-2, 1:-1])/(2.*Grid.dz)
      fy[1:-1, 1:-1] = (phi[1:-1, 2:] - phi[1:-1, :-2])/(2.*Grid.dy)
     
      #print np.mean(par.dt*fz*var.rho), np.mean(var.momentumZ)

      momentumGravSourceZ = 0.5*(var.rho+var.lastrho)*fz
      momentumGravSourceY = 0.5*(var.rho+var.lastrho)*fy
      energyGravSource = var.momentumZ*fz + var.momentumY*fy

   elif par.GravityMode == 'PoissonImplicit':
      phi = np.zeros(Grid.z.shape)
      fz = np.zeros(Grid.z.shape)
      fy = np.zeros(Grid.z.shape)

      nnz, col, offset = BuildSystemMatrix(-2*(1./(Grid.dz*Grid.dz) + 1./(Grid.dy*Grid.dy)), 1./(Grid.dy*Grid.dy), 1./(Grid.dy*Grid.dy), 1./(Grid.dz*Grid.dz), 1./(Grid.dz*Grid.dz) )
      N = Grid.z.shape[0]*Grid.z.shape[1]   #Repeated operation
      SystemMatrix = scipy.sparse.csr_matrix( (nnz, col, offset), shape=(N, N)  )
      
      #var.rho = 0.1*np.ones(Grid.z.shape)*np.exp( -( (Grid.z-0.4)*(Grid.z-0.4)+(Grid.y-0.4)*(Grid.y-0.4) )/0.005 )
      #var.rho += 0.1*np.exp( -( (Grid.z-0.7)*(Grid.z-0.7)+(Grid.y-0.7)*(Grid.y-0.7) )/0.005 )
      
      sol = scipy.sparse.linalg.spsolve(SystemMatrix, var.rho.flatten() )
      phi = sol.reshape(Grid.z.shape)
      
      fy[1:-1, 1:-1] = -(phi[2:, 1:-1] - phi[:-2, 1:-1])/(2.*Grid.dz)
      fz[1:-1, 1:-1] = -(phi[1:-1, 2:] - phi[1:-1, :-2])/(2.*Grid.dy)

      
      """
      import matplotlib.pyplot as plt
      plt.pcolor(Grid.z[1:-1,1:-1], Grid.y[1:-1,1:-1], phi[2:, 1:-1]-phi[:-2, 1:-1])
      plt.colorbar()
      
      plt.contour(Grid.z, Grid.y, var.rho)
      
      plt.quiver(Grid.z, Grid.y, fz, fy)
      plt.show()
      """
      momentumGravSourceZ = (var.rho)*fz
      momentumGravSourceY = (var.rho)*fy
      energyGravSource = var.momentumZ*fz + var.momentumY*fy

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

   return momentumDampingSource, energyDampingSource

def computeTemperatureDiffusion():
   if par.SpitzerDiffusion:
      var.kappa = par.ct*var.T**(5./2.)
   
   if par.dim == 1:
      Dkappa = (var.kappa[2:]-var.kappa[:-2])/4.
      EnergyDiff = (var.kappa[1:-1]+Dkappa)*var.T[2:] - 2.*var.kappa[1:-1]*var.T[1:-1] + (var.kappa[1:-1]-Dkappa)*var.T[:-2]

      return par.DiffusionPercent*EnergyDiff/(Grid.dz*Grid.dz)
   
   elif par.dim==2:
      DkappaZ = (var.kappa[1:-1, 2:]-var.kappa[1:-1, :-2])/4.
      DkappaY = (var.kappa[2:, 1:-1]-var.kappa[:-2, 1:-1])/4.
   
      dz2 = Grid.dz*Grid.dz
      dy2 = Grid.dy*Grid.dy

      EnergyDiffZ = (var.kappa[1:-1,1:-1]+DkappaZ)*var.T[1:-1,2:] - 2.*var.kappa[1:-1,1:-1]*var.T[1:-1,1:-1] + (var.kappa[1:-1,1:-1]-DkappaZ)*var.T[1:-1,:-2]
      EnergyDiffY = (var.kappa[1:-1,1:-1]+DkappaY)*var.T[2:,1:-1] - 2.*var.kappa[1:-1,1:-1]*var.T[1:-1,1:-1] + (var.kappa[1:-1,1:-1]-DkappaY)*var.T[:-2,1:-1]
      
      return EnergyDiffZ/dz2 + EnergyDiffY/dy2

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
   if par.dim == 1:
      e = par.cv*var.T   #Internal energy
      Ek = var.energy-e*var.rho
      Dkappa = (var.kappa[2:]-var.kappa[:-2])/4.
   
      dz2 = Grid.dz*Grid.dz
      A = par.DiffusionPercent * (var.kappa[1:-1]+Dkappa)/dz2
      B = -2.*par.DiffusionPercent * var.kappa[1:-1]/dz2
      C = par.DiffusionPercent * (var.kappa[1:-1]-Dkappa)/dz2
   
      var.diag[1:-1] = par.dt*B/(var.rho[1:-1]*par.cv) - 1.  
      var.lower[:-2] = par.dt*C/(var.rho[:-2]*par.cv)   
      var.upper[2:] = par.dt*A/(var.rho[2:]*par.cv)
   
      var.rhs[1:-1] = -e[1:-1]*var.rho[1:-1]

      if sets.BoundaryConditionL!=None: 
         sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ,), var.energy)
      if sets.BoundaryConditionR!=None: 
         sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ,), var.energy)
   
   
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
         var.energy = sol_e + Ek #0.5*var.momentumZ*var.momentumZ/var.rho
      elif par.dim == 2:
         var.energy = sol_e + 0.5*( var.momentumZ*var.momentumZ + var.momentumY*var.momentumY )/var.rho
   
      var.T = sol_e/(var.rho*par.cv)

   elif par.dim == 2:
      e = par.cv*var.T
      Ek = var.energy-var.rho*e

      Dkappa = 0.5*(var.kappa[1:-1,2:] - var.kappa[1:-1,:-2])/Grid.dz + 0.5*(var.kappa[2:,1:-1] - var.kappa[:-2,1:-1])/Grid.dy
              
      dz2 = Grid.dz*Grid.dz
      dy2 = Grid.dy*Grid.dy
      
      
      var.rhs = np.zeros_like(Grid.z)
      var.rhs[1:-1,1:-1] = -e[1:-1,1:-1]*var.rho[1:-1,1:-1]
      var.rhs = var.rhs.flatten()
      
      L = np.zeros_like(Grid.z)
      R = np.zeros_like(Grid.z)
      C = np.zeros_like(Grid.z)
      T = np.zeros_like(Grid.z)
      B = np.zeros_like(Grid.z)

      C[1:-1, 1:-1] = par.dt * var.kappa[1:-1,1:-1] * (-2./dy2  - 2./dz2) / (var.rho[1:-1,1:-1]*par.cv) - 1.
      R[1:-1, 1:-1] = par.dt * ( 0.5*Dkappa/Grid.dz + var.kappa[1:-1,1:-1]/dz2 ) / (var.rho[1:-1,2:]*par.cv)
      L[1:-1, 1:-1] = par.dt * (-0.5*Dkappa/Grid.dz + var.kappa[1:-1,1:-1]/dz2 ) / (var.rho[1:-1,:-2]*par.cv)
      T[1:-1, 1:-1] = par.dt * ( 0.5*Dkappa/Grid.dy + var.kappa[1:-1,1:-1]/dy2 ) / (var.rho[2:,1:-1]*par.cv)
      B[1:-1, 1:-1] = par.dt * (-0.5*Dkappa/Grid.dy + var.kappa[1:-1,1:-1]/dy2 ) / (var.rho[:-2,1:-1]*par.cv)
   
      #var.diag[1:-1] = par.dt*B/(var.rho[1:-1]*par.cv) - 1.  
      #var.lower[:-2] = par.dt*C/(var.rho[:-2]*par.cv)   
      #var.upper[2:] = par.dt*A/(var.rho[2:]*par.cv)

      nnz, col, offset = BuildSystemMatrix(C, T, B, L, R) 
      N = Grid.z.shape[0]*Grid.z.shape[1]   
      
      #print nnz, col, offset, nnz.shape, offset.shape
      #var.rho = 0.1*np.ones(Grid.z.shape)*np.exp( -( (Grid.z-0.4)*(Grid.z-0.4)+(Grid.y-0.4)*(Grid.y-0.4) )/0.005 )
      #var.rho += 0.1*np.exp( -( (Grid.z-0.7)*(Grid.z-0.7)+(Grid.y-0.7)*(Grid.y-0.7) )/0.005 )
      
      
      
      #print rhs.flatten().shape, N
      """
      if sets.BoundaryConditionL!=None: 
         sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ,), var.energy)
      if sets.BoundaryConditionR!=None: 
         sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ,), var.energy)
      """
      
      SystemMatrix = scipy.sparse.csr_matrix( (nnz, col, offset), shape=(N, N)  )
      
      sol_e, info = scipy.sparse.linalg.gmres(SystemMatrix, var.rhs.flatten(), x0=var.T.flatten() )
      sol_e = sol_e.reshape(Grid.z.shape)
      if info!=0:
         print '### GMRES DID NOT CONVERGE'
      
   
      #var.rhs[1:-1] = -e[1:-1]*var.rho[1:-1]


      
      #print e[:6]*var.rho[:6]
      #sol_e = linalg.solve_banded( (1,1), [var.upper, var.diag, var.lower], var.rhs )
      #print sol_e[:6]
      #print ''
      #print var.T-sol_e/(par.cv*var.rho)
      
      var.T = sol_e/(var.rho*par.cv)
      var.energy = Ek + sol_e
   
            




def computeRadiativeLosses():
   logT = np.log10(var.T)

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
      if par.staggered:
         var.momentumZ[1:-1, :-1] += par.dt*momentumGZ
         var.momentumY += par.dt*momentumGY
         var.energy[1:-1, 1:-1] += par.dt*energyG        
      else: 
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
         if par.dim == 1:
            var.energy[1:-1] += par.dt*ThermalDiff
         elif par.dim == 2:
            var.energy[1:-1,1:-1] += par.dt*ThermalDiff
            
         ChangeOfVar.ConvertToPrim()   

   if par.MomentumDamping:
      momentumDamping, energyDamping = computeMomentumDamping()
      var.momentum += momentumDamping
      var.energy += energyDamping
      ChangeOfVar.ConvertToPrim()
     
     
     
     
     
     
     
