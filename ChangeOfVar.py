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
   var.v = var.momentum/var.rho
   var.P = (var.energy - 0.5*var.momentum*var.momentum/var.rho)*(par.gamma-1.)
   var.T = var.P*par.mu/(var.rho*par.R)
   
   #e = var.energy/var.rho - 0.5*var.v*var.v              #rho*e = E - 0.5*rho*v*v

   #var.T = e/par.cv          
   #var.P = var.rho*(par.gamma-1.)*e
   
