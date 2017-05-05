import numpy as np

import Variables as var
import matplotlib.pyplot as plt

# ------------------ PLOT ----------------------

#fig = plt.figure()
#ax = fig.add_subplot(111)
#line1, = ax.plot(z, rho)


f, axs = plt.subplots(3,1, sharex=True)
f.set_size_inches(18.5, 10.5)
axs[0].set_title(r"$\rho$")
rho_line, = axs[0].plot(var.z, var.rho)
axs[1].set_title("v")
v_line, = axs[1].plot(var.z, var.v)
axs[2].set_title("P")
P_line, = axs[2].plot(var.z, var.P)

FreeFall_line, = axs[1].plot([],[], "r--")
SoundSpeed_line, = axs[1].plot([],[], "r--")

if var.IsThereGravity:
   FreeFall_line.set_xdata([var.z[0], var.z[-1]])
if var.SoundSpeedLine:
   SoundSpeed_line.set_ydata([0.,2.])

axs[0].set_ylim(0.995,1.005)
axs[1].set_ylim(-.005, .005)
axs[2].set_ylim(0.995, 1.005)
axs[0].set_xlim(var.z0,var.zf)

def Plot():
   rho_line.set_ydata(var.rho)
   v_line.set_ydata(var.v)
   P_line.set_ydata(var.P)

   if var.IsThereGravity:
      vff = var.g*var.tt
      FreeFall_line.set_ydata([vff, vff])

   if var.SoundSpeedLine:
      c_s = np.sqrt(var.gamma)   #rho=p0=1
      x = c_s*var.tt 
      SoundSpeed_line.set_xdata([x, x])

   plt.savefig('RESULTS/%.5f.png'%var.tt, bbox_inches='tight')


