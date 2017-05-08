import Parameters as par
import Variables as var

print "Loading ChangeOfpar.."

def ConvertToPrim():
   var.v = var.momentum/var.rho
   e = var.energy/var.rho - 0.5*var.v*var.v              #rho*e = E - 0.5*rho*v*v
   
   var.T = e/par.Cv          
   var.P = var.rho*(par.gamma-1.)*e
