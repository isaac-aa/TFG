#--------------------------------------------
#		ChangeOfVar.py
# This module changes from Conservative (rho,
# momentum, energy) to Primitive variables
# (rho, velocity, pressure) 
#
#--------------------------------------------


import Parameters as par
import Variables as var
import Grid
print "Loading ChangeOfpar.."


def ConvertToPrim():
   if not par.staggered:
      var.vZ = var.momentumZ/var.rho
      var.vY = var.momentumY/var.rho
      var.P = (var.energy - 0.5*(var.momentumZ*var.momentumZ + var.momentumY*var.momentumY)/var.rho)*(par.gamma-1.)
      var.T = var.P*par.mu/(var.rho*par.R)
      
   else:  
      if par.dim==1:
         var.vZ[:-1] = 2.*var.momentumZ[:-1]/(var.rho[:-1] + var.rho[1:])
         #var.vY[:-1] = 2.*var.momentumY[:-1]/(var.rho[:-1] + var.rho[1:])
         momZ = 0.5*(var.momentumZ[:-2]+var.momentumZ[1:-1])
         var.P[1:-1] = (var.energy[1:-1] - 0.5*(momZ*momZ)/var.rho[1:-1])*(par.gamma-1.)
         var.P[0] = (var.energy[0] - 0.5*(var.momentumZ[0]*var.momentumZ[0])/var.rho[0])*(par.gamma-1.)
         var.P[-1] = (var.energy[-1] - 0.5*(var.momentumZ[-2]*var.momentumZ[-2])/var.rho[-1])*(par.gamma-1.)
         var.T = var.P*par.mu/(var.rho*par.R)
         
      elif par.dim==2:
         var.vZ[1:-1,:-1] = 2.*var.momentumZ[1:-1,:-1]/(var.rho[1:-1,:-1] + var.rho[1:-1,1:])
         var.vY[:-1,1:-1] = 2.*var.momentumY[:-1,1:-1]/(var.rho[:-1,1:-1] + var.rho[1:,1:-1])
         momZ = 0.5*(var.momentumZ[1:-1,:-2]+var.momentumZ[1:-1,1:-1])
         momY = 0.5*(var.momentumY[:-2,1:-1]+var.momentumY[1:-1,1:-1])
         var.P[1:-1,1:-1] = (var.energy[1:-1,1:-1] - 0.5*(momZ*momZ+momY*momY)/var.rho[1:-1,1:-1])*(par.gamma-1.)
         var.P[1:-1,   0] = (var.energy[1:-1,0] - 0.5*(momZ[:,0]*momZ[:,0]+var.momentumY[1:-1, 0]*var.momentumY[1:-1, 0])/var.rho[1:-1,0])*(par.gamma-1.) # L
         var.P[1:-1,  -1] = (var.energy[1:-1,-1] - 0.5*(momZ[:,-1]*momZ[:,-1]+var.momentumY[1:-1,-2]*var.momentumY[1:-1,-2])/var.rho[1:-1,-1])*(par.gamma-1.) # R
         var.P[0,   1:-1] = (var.energy[0,1:-1] - 0.5*(var.momentumZ[0,1:-1]*var.momentumZ[0,1:-1]+momY[0,:]*momY[0,:])/var.rho[0,1:-1])*(par.gamma-1.) # B
         var.P[-1,  1:-1] = (var.energy[-1,1:-1] - 0.5*(var.momentumZ[-2,1:-1]*var.momentumZ[-2,1:-1]+momY[-1,:]*momY[-1,:])/var.rho[-1,1:-1])*(par.gamma-1.) # T
         if par.MHD==True:
            var.P += (- 0.5*(var.Bx*var.Bx + var.By*var.By + var.Bz*var.Bz) - 0.5*var.momentumX*var.momentumX/var.rho )* (par.gamma-1.)
            var.vX = var.momentumX/var.rho
         
         var.T = var.P*par.mu/(var.rho*par.R)        


def ConvertToPrimBoundaries():
   if not par.staggered:
      var.vZ[0] = var.momentumZ/var.rho
      var.vY[0] = var.momentumY/var.rho
      var.P[0] = (var.energy - 0.5*(var.momentumZ*var.momentumZ + var.momentumY*var.momentumY)/var.rho)*(par.gamma-1.)
      var.T[0] = var.P*par.mu/(var.rho*par.R)
      
      var.vZ[-1] = var.momentumZ/var.rho
      var.vY[-1] = var.momentumY/var.rho
      var.P[-1] = (var.energy - 0.5*(var.momentumZ*var.momentumZ + var.momentumY*var.momentumY)/var.rho)*(par.gamma-1.)
      var.T[-1] = var.P*par.mu/(var.rho*par.R)
      
   else:  
      if par.dim==1:
         var.vZ[-2] = 2.*var.momentumZ[:-1]/(var.rho[:-1] + var.rho[1:])
         #var.vY[:-1] = 2.*var.momentumY[:-1]/(var.rho[:-1] + var.rho[1:])
         var.P[0] = (var.energy[0] - 0.5*(var.momentumZ[0]*var.momentumZ[0])/var.rho[0])*(par.gamma-1.)
         var.P[-1] = (var.energy[-1] - 0.5*(var.momentumZ[-2]*var.momentumZ[-2])/var.rho[-1])*(par.gamma-1.)
         var.T[-1] = var.P[-1]*par.mu/(var.rho[-1]*par.R)

         var.vZ[0] = 2.*var.momentumZ[0]/(var.rho[0] + var.rho[1])
         var.T[0] = var.P[0]*par.mu/(var.rho[0]*par.R)

      
      elif par.dim==2:
         
         var.vZ[1:-1,   0] = 2.*var.momentumZ[1:-1,   0]/(var.rho[1:-1,   0] + var.rho[1:-1,   1])
         var.vZ[1:-1,  -2] = 2.*var.momentumZ[1:-1,  -2]/(var.rho[1:-1,  -2] + var.rho[1:-1,  -1])
         var.vZ[0,   1:-1] = 2.*var.momentumZ[0,   1:-1]/(var.rho[0,   :-2] + var.rho[0,   1:-1])
         var.vZ[-1,  1:-1] = 2.*var.momentumZ[-1,  1:-1]/(var.rho[-1,  1:-1] + var.rho[-1,  1:-1])
         

         var.vY[1:-1,   0] = 2.*var.momentumZ[1:-1,   0]/(var.rho[:-2,   0] + var.rho[1:-1,   0])
         var.vY[1:-1,  -1] = 2.*var.momentumZ[1:-1,  -1]/(var.rho[:-2,  -1] + var.rho[1:-1,  -1])
         var.vY[0,   1:-1] = 2.*var.momentumZ[0,   1:-1]/(var.rho[0,   1:-1] + var.rho[1,   1:-1])
         var.vY[-2,  1:-1] = 2.*var.momentumZ[-2,  1:-1]/(var.rho[-1,  1:-1] + var.rho[-2,  1:-1])
         
         """
         momZ0 = 0.5*(var.momentumZ[1:-1,-3]+var.momentumZ[1:-1,-2])  #Old but good-looking version
         momY0 = 0.5*(var.momentumY[-3,1:-1]+var.momentumY[-2,1:-1])
         momZ1 = var.momentumZ[1:-1,0]  #0.5*(var.momentumZ[1:-1, 0]+var.momentumZ[1:-1, 1])
         momY1 = 0.5*(var.momentumY[0 ,1:-1]+var.momentumY[1 ,1:-1])
 
         var.P[1:-1,   0] = (var.energy[1:-1,0] - 0.5*(momZ1*momZ1+var.momentumY[1:-1, 0]*var.momentumY[1:-1, 0])/var.rho[1:-1,0])*(par.gamma-1.)# L
         var.P[1:-1,  -1] = (var.energy[1:-1,-1] - 0.5*(momZ0*momZ0+var.momentumY[1:-1,-2]*var.momentumY[1:-1,-2])/var.rho[1:-1,-1])*(par.gamma-1.) # R
         var.P[0,   1:-1] = (var.energy[0,1:-1] - 0.5*(var.momentumZ[0,1:-1]*var.momentumZ[0,1:-1]+momY1*momY1)/var.rho[0,1:-1])*(par.gamma-1.) # B
         var.P[-1,  1:-1] = (var.energy[-1,1:-1] - 0.5*(var.momentumZ[-2,1:-1]*var.momentumZ[-2,1:-1]+momY0*momY0)/var.rho[-1,1:-1])*(par.gamma-1.) # T
         """
         
         momYleft = 0.5*(var.momentumY[:-2,0] + var.momentumY[1:-1,0])
         momYright = 0.5*(var.momentumY[:-2,-1] + var.momentumY[1:-1,-1])
         momZtop = 0.5*(var.momentumZ[-1,:-2] + var.momentumZ[-1,1:-1])
         momZbot = 0.5*(var.momentumZ[0,:-2] + var.momentumZ[0,1:-1])
         
         # TEMPORARY CHANGE
         momZright = 0.5*(Grid.extraZ[1:-1] + var.momentumZ[1:-1, 0])
         momZleft = 0.5*(var.momentumZ[1:-1,-2] + var.momentumZ[1:-1, -1])
         momYtop = 0.5*(var.momentumY[-1,1:-1] + var.momentumY[-2, 1:-1])
         momYbot = 0.5*(Grid.extraY[1:-1] + var.momentumY[0, 1:-1])

         var.P[1:-1,   0] = (var.energy[1:-1,0] - 0.5*(momZleft*momZleft + momYleft*momYleft)/var.rho[1:-1,0])*(par.gamma-1.)# L
         var.P[1:-1,  -1] = (var.energy[1:-1,-1] - 0.5*(momZright*momZright + momYright*momYright)/var.rho[1:-1,-1])*(par.gamma-1.) # R
         var.P[0,   1:-1] = (var.energy[0,1:-1] - 0.5*(momZbot*momZbot + momYbot*momYbot)/var.rho[0,1:-1])*(par.gamma-1.) # B
         var.P[-1,  1:-1] = (var.energy[-1,1:-1] - 0.5*(momZtop*momZtop + momYtop*momYtop)/var.rho[-1,1:-1])*(par.gamma-1.) # T  

         
         
         #var.P[1:-1,   0] = (var.energy[1:-1,0] - 0.5*(var.momentumZ[1:-1,0]*var.momentumZ[1:-1,0] + momYleft*momYleft)/var.rho[1:-1,0])*(par.gamma-1.)# L
         #var.P[1:-1,  -1] = (var.energy[1:-1,-1] - 0.5*(var.momentumZ[1:-1,-2]*var.momentumZ[1:-1,-2] + momYright*momYright)/var.rho[1:-1,-1])*(par.gamma-1.) # R
         #var.P[0,   1:-1] = (var.energy[0,1:-1] - 0.5*(momZbot*momZbot + var.momentumY[0,1:-1]*var.momentumY[0,1:-1])/var.rho[0,1:-1])*(par.gamma-1.) # B
         #var.P[-1,  1:-1] = (var.energy[-1,1:-1] - 0.5*(momZtop*momZtop + var.momentumY[-2,1:-1]*var.momentumY[-2,1:-1])/var.rho[-1,1:-1])*(par.gamma-1.) # T  

         
         if par.MHD:
            var.P[1:-1,   0] += -(0.5*(var.Bx[1:-1,   0]*var.Bx[1:-1,   0] + var.By[1:-1,   0]*var.By[1:-1,   0] + var.Bz[1:-1,   0]*var.Bz[1:-1,   0]) + 0.5*(var.momentumX[1:-1,   0]*var.momentumX[1:-1,   0]/var.rho[1:-1,   0]))*(par.gamma-1.)
            var.P[1:-1,  -1] += -(0.5*(var.Bx[1:-1,  -1]*var.Bx[1:-1,  -1] + var.By[1:-1,  -1]*var.By[1:-1,  -1] + var.Bz[1:-1,  -1]*var.Bz[1:-1,  -1]) + 0.5*(var.momentumX[1:-1,  -1]*var.momentumX[1:-1,  -1]/var.rho[1:-1,  -1]))*(par.gamma-1.)
            var.P[0,   1:-1] += -(0.5*(var.Bx[0,   1:-1]*var.Bx[0,   1:-1] + var.By[0,   1:-1]*var.By[0,   1:-1] + var.Bz[0,   1:-1]*var.Bz[0,   1:-1]) + 0.5*(var.momentumX[0,   1:-1]*var.momentumX[0,   1:-1]/var.rho[0,   1:-1]))*(par.gamma-1.)
            var.P[-1,  1:-1] += -(0.5*(var.Bx[-1,  1:-1]*var.Bx[-1,  1:-1] + var.By[-1,  1:-1]*var.By[-1,  1:-1] + var.Bz[-1,  1:-1]*var.Bz[-1,  1:-1]) + 0.5*(var.momentumX[-1,  1:-1]*var.momentumX[-1,  1:-1]/var.rho[-1,  1:-1]))*(par.gamma-1.)
            
            var.vX[1:-1,   0] = var.momentumX[1:-1,   0]/var.rho[1:-1,   0]
            var.vX[1:-1,  -1] = var.momentumX[1:-1,  -1]/var.rho[1:-1,  -1]
            var.vX[0,   1:-1] = var.momentumX[0,   1:-1]/var.rho[0,   1:-1]
            var.vX[-1,  1:-1] = var.momentumX[-1,  1:-1]/var.rho[-1,  1:-1]

         var.T[1:-1,   0] = var.P[1:-1,   0]*par.mu/(var.rho[1:-1,   0]*par.R)        
         var.T[1:-1,  -1] = var.P[1:-1,  -1]*par.mu/(var.rho[1:-1,  -1]*par.R)        
         var.T[0,   1:-1] = var.P[0,   1:-1]*par.mu/(var.rho[0,   1:-1]*par.R)        
         var.T[-1,  1:-1] = var.P[-1,  1:-1]*par.mu/(var.rho[-1,  1:-1]*par.R)        

