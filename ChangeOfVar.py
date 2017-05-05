import Variables as var

def ConvertToPrim(rho, momentum, energy):
   v = momentum/rho
   e = energy/rho - 0.5*v*v              #rho*e = E - 0.5*rho*v*v
   
   T = e/var.Cv          
   P = rho*(var.gamma-1.)*e
   return v,T,P
