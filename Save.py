#--------------------------------------------
#		Save.py
# This module saves the current state as a plot,
# or a ASCII file (TODO). It can add different
# plots, that can be set at Parameters.py (eg.
# an analytic solution)
#
#--------------------------------------------

import numpy as np

import Grid
import Parameters as par
import Variables as var
import Settings as sets
import Analytic
import matplotlib.pyplot as plt

# ------------------ PLOT ----------------------

#fig = plt.figure()
#ax = fig.add_subplot(111)
#line1, = ax.plot(z, rho)


f, axs = plt.subplots(4,1, sharex=True)
f.set_size_inches(18.5, 10.5)
axs[0].set_title(r"$\rho$")
rho_line, = axs[0].plot(Grid.z, var.rho)
axs[1].set_title("v")
v_line, = axs[1].plot(Grid.z, var.v)
axs[2].set_title("P")
P_line, = axs[2].plot(Grid.z, var.P)
axs[3].set_title("T")
T_line, = axs[3].plot(Grid.z, var.T)


FreeFall_line, = axs[1].plot([],[], "r--")

SoundSpeed_line, = axs[1].plot([],[], "r--")
MaxAmp_line, = axs[1].plot([],[], "k--")

rhoAna_line, = axs[0].plot([],[], "k--")
vAna_line, = axs[1].plot([],[], "k--")
PAna_line, = axs[2].plot([],[], "k--")


if par.IsThereGravity:
   FreeFall_line.set_xdata([Grid.z[0], Grid.z[-1]])
if par.SoundSpeedLine:
   SoundSpeed_line.set_ydata([0.,2.])
   MaxAmp_line.set_xdata([0.,1.])
if par.SoundSpeedAnalytic or par.IsothermalAnalytic:
   rhoAna_line.set_xdata(Grid.z)
   vAna_line.set_xdata(Grid.z)
   PAna_line.set_xdata(Grid.z)

if par.rhoAxis != []:
  axs[0].set_ylim(par.rhoAxis)
if par.vAxis != []:
  axs[1].set_ylim(par.vAxis)
if par.PAxis != []:
  axs[2].set_ylim(par.PAxis)
if par.TAxis != []:
  axs[3].set_ylim(par.TAxis)
  
for i in range(4):
   if par.logScale[i]==True:
      axs[i].semilogy()

axs[0].set_xlim(Grid.z[0],Grid.z[-1])

def Plot():
   rho_line.set_ydata(var.rho)
   v_line.set_ydata(var.v)
   P_line.set_ydata(var.P)
   T_line.set_ydata(var.T)

   if par.IsThereGravity:
      vff = par.g*par.tt
      FreeFall_line.set_ydata([vff, vff])

   if par.SoundSpeedLine:
      c_s = np.sqrt(par.gamma)   #rho=p0=1
      x = c_s*par.tt 
      SoundSpeed_line.set_xdata([x, x])
      MaxAmp = np.max(var.v)
      MaxAmp_line.set_ydata([MaxAmp, MaxAmp])
      
      
   # This can be further improved...
   if par.SoundSpeedAnalytic:
      rhoAna, vAna, PAna = Analytic.SoundWaves(sets.argsIC, par.tt)
      rhoAna_line.set_ydata(rhoAna)
      vAna_line.set_ydata(vAna)
      PAna_line.set_ydata(PAna)
   
   if par.IsothermalAnalytic:
      rhoAna, vAna, PAna = Analytic.Isothermal(sets.argsIC, par.tt)
      rhoAna_line.set_ydata(rhoAna)
      vAna_line.set_ydata(vAna)
      PAna_line.set_ydata(PAna)
   
   
   plt.savefig('RESULTS/%.5f.png'%par.tt, bbox_inches='tight')


