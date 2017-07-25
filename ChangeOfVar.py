#--------------------------------------------
#		ChangeOfVar.py
# This module changes from Conservative (rho,
# momentum, energy) to Primitives variables
# (rho, velocity, pressure, temperature) and
# viceversa (TODO)
#
#--------------------------------------------


import Parameters as par
import Variables as var

print "Loading ChangeOfpar.."

def ConvertToPrim():
   var.vZ = var.momentumZ/var.rho
   var.vY = var.momentumY/var.rho
   e = var.energy/var.rho - 0.5*(var.vZ*var.vZ + var.vY*var.vY)              #rho*e = E - 0.5*rho*v*v

   var.T = e/par.cv          
   var.P = var.rho*(par.gamma-1.)*e
   
