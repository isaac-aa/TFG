#--------------------------------------------
#		ChangeOfVar.py
# This module changes from Conservative (rho,
# momentum, energy) to Primitive variables
# (rho, velocity, pressure) 
#
#--------------------------------------------


import Parameters as par
import Variables as var

print "Loading ChangeOfpar.."

def ConvertToPrim():
   if not par.staggered:
      var.vZ = var.momentumZ/var.rho
      var.vY = var.momentumY/var.rho
      var.P = (var.energy - 0.5*(var.momentumZ*var.momentumZ + var.momentumY*var.momentumY)/var.rho)*(par.gamma-1.)
      var.T = var.P*par.mu/(var.rho*par.R)
      
   else:  #Only for 1D
      var.vZ[:-1] = 2.*var.momentumZ[:-1]/(var.rho[:-1] + var.rho[1:])
      var.vY[:-1] = 2.*var.momentumY[:-1]/(var.rho[:-1] + var.rho[1:])
      momZ = 0.5*(var.momentumZ[:-2]+var.momentumZ[1:-1])
      var.P[1:-1] = (var.energy[1:-1] - 0.5*(momZ*momZ)/var.rho[1:-1])*(par.gamma-1.)
      var.P[0] = (var.energy[0] - 0.5*(var.momentumZ[0]*var.momentumZ[0])/var.rho[0])*(par.gamma-1.)
      var.P[-1] = (var.energy[-1] - 0.5*(var.momentumZ[-2]*var.momentumZ[-2])/var.rho[-1])*(par.gamma-1.)
      var.T = var.P*par.mu/(var.rho*par.R)
   
