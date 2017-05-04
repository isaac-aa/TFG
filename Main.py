import numpy as np
import time
import matplotlib.pyplot as plt
from Variables import * 
import SourceTerm
from Settings import *

# ------------------ PLOT ----------------------

#fig = plt.figure()
#ax = fig.add_subplot(111)
#line1, = ax.plot(z, rho)

f, axs = plt.subplots(3,1, sharex=True)
f.set_size_inches(18.5, 10.5)
axs[0].set_title(r"$\rho$")
rho_line, = axs[0].plot(z, rho)
axs[1].set_title("v")
v_line, = axs[1].plot(z, momentum/rho)
axs[2].set_title("P")
P_line, = axs[2].plot(z, T*rho)

FreeFall_line, = axs[1].plot([],[], "r--")

if IsThereGravity:
   FreeFall_line.set_xdata([z[0], z[-1]])

axs[0].set_ylim(0.95,1.05)
axs[1].set_ylim(-.05, .05)
axs[2].set_ylim(0.95, 1.05)


def Plot(t):
   rho_line.set_ydata(rho)
   v_line.set_ydata(v)
   P_line.set_ydata(P)

   if IsThereGravity:
      vff = g*t
      FreeFall_line.set_ydata([vff, vff])

   plt.savefig('RESULTS/%.5f.png'%t, bbox_inches='tight')


# ------------------ TIME STEP -----------------





while (it<max_it and tt<tf):

   # Boundary conditions
   BoundaryConditionL(rho, momentum, energy, argsL)
   BoundaryConditionR(rho, momentum, energy, argsR)

   # Source computation
   if IsComputingSource:
      momentumSource, energySource = SourceTerm.ComputeSource(rho, momentum)

   # Time step
   rho, momentun, energy = FluxScheme(dt, dz, rho, momentum, momentumSource, energy, energySource) 

   # Compute change of variables
   v = momentum/rho
   
   e = energy/rho - 0.5*v*v              #rho*e = E - 0.5*rho*v*v
   
   T = e/Cv          
   P = rho*(gamma-1.)*e



   if it%500 == 0.:
     Plot(tt)
  
   
   cfl = dt*np.max(abs(momentum/rho))/np.min(dz)
   if np.max(abs(momentum))!=0:
     dt = cfl_set*np.min(dz)/np.max(abs(momentum/rho))
 
   if dt>=dt_max:
      dt = dt_max

   print it, tt, dt,  cfl
   tt = tt+dt
   it = it+1
   #time.sleep(1)






