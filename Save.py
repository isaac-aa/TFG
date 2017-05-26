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
import Characteristics

# ------------------ PLOT ----------------------

#fig = plt.figure()
#ax = fig.add_subplot(111)
#line1, = ax.plot(z, rho)
print 'Loading Save..'

plotCounter = 0

f, axs = plt.subplots(4+par.PlotCharacteristics,1, sharex=True)
f.set_size_inches(18.5, 10.5)
axs[0].set_title(r"$\rho$")
rho_line, = axs[0].plot(Grid.z, var.rho)
axs[1].set_title("v")
v_line, = axs[1].plot(Grid.z, var.v)
axs[2].set_title("P")
P_line, = axs[2].plot(Grid.z, var.P)
axs[3].set_title("T")
T_line, = axs[3].plot(Grid.z, var.T)

if par.PlotCharacteristics:
  axs[4].set_title('Characteristic lenght/time')
  ax2 = axs[4].twinx()
  if par.ThermalDiffusion:
    Thermal_tau, = axs[4].plot([],[], 'k')
    Thermal_tau.set_xdata(Grid.z)
    Thermal_L, = ax2.plot([],[], 'k--')
    Thermal_L.set_xdata(Grid.z) 
  if par.RadiativeLoss:
    Losses_tau, = axs[4].plot([],[], 'g')
    Losses_tau.set_xdata(Grid.z)
  Dynamic_tau, = axs[4].plot([],[], 'b')
  Dynamic_tau.set_xdata(Grid.z)
  axs[4].semilogy()
  ax2.semilogy()

FreeFall_line, = axs[1].plot([],[], "r--")

SoundSpeed_line, = axs[1].plot([],[], "r--")
MaxAmp_line, = axs[1].plot([],[], "k--")

rhoAna_line, = axs[0].plot([],[], "k--")
vAna_line, = axs[1].plot([],[], "k--")
PAna_line, = axs[2].plot([],[], "k--")
TAna_line, = axs[3].plot([],[], "k--")

SoundSpeedProf_line, = axs[1].plot([],[], 'g--')

if par.IsThereGravity:
   FreeFall_line.set_xdata([Grid.z[0], Grid.z[-1]])

if par.SoundSpeedLine:
   SoundSpeed_line.set_ydata([0.,2.])
   MaxAmp_line.set_xdata([0.,1.])

if par.SoundSpeedProfile:
   SoundSpeedProf_line.set_xdata(Grid.z)
   
if par.SoundSpeedAnalytic or par.IsothermalAnalytic or par.ThermalAnalytic:
   rhoAna_line.set_xdata(Grid.z)
   vAna_line.set_xdata(Grid.z)
   PAna_line.set_xdata(Grid.z)
   TAna_line.set_xdata(Grid.z)


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


if par.PlotFile:
   z_dat, rho_dat, p_dat, T_dat = np.loadtxt(par.FileToPlot, unpack=True)
   axs[0].plot(z_dat, rho_dat, 'k--')
   axs[2].plot(z_dat, p_dat, 'k--')
   axs[3].plot(z_dat, T_dat, 'k--')


def Plot():
   global plotCounter

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
      
   if par.SoundSpeedProfile:
      c_s = np.sqrt(par.gamma*var.P/var.rho)
      SoundSpeedProf_line.set_ydata(c_s)
      
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
      
   if par.ThermalAnalytic:
      rhoAna, vAna, PAna, TAna = Analytic.GaussianThermal(sets.argsIC, par.tt)
      rhoAna_line.set_ydata(rhoAna)
      vAna_line.set_ydata(vAna)
      PAna_line.set_ydata(PAna)
      TAna_line.set_ydata(TAna)
      
   if par.PlotCharacteristics:
     if par.ThermalDiffusion:
       tau, L = Characteristics.Thermal()
       Thermal_tau.set_ydata(tau)
       Thermal_L.set_ydata(L) 
     if par.RadiativeLoss:
       tau = Characteristics.Radiative()
       Losses_tau.set_ydata(tau)
    
     tau = Characteristics.DensityChanges()
     Dynamic_tau.set_ydata(tau) 
     print 'Min characteristic time: %.2f'%np.min(np.abs(tau))
     #Re-scale axis
     axs[4].relim()
     axs[4].autoscale_view()  
     ax2.relim()
     ax2.autoscale_view()  
   #print '----' 
   #print var.T[0], var.T[1], (var.T[0]+var.T[1])/2.
   #print var.T[-2], var.T[-1], (var.T[-1]+var.T[-2])/2.
   #print '----' 
   if par.SaveToFile and plotCounter%par.SaveToFileRatio==0 :
      dataToSave = np.array([Grid.z, var.rho, var.v, var.T])
      if par.PlotCharacteristics:  
        if par.ThermalDiffusion:
           tau_T = Thermal_tau.get_ydata()
           L_T = Thermal_L.get_ydata()
           dataToSave = np.append(dataToSave, [tau_T, L_T], axis=0)
        if par.RadiativeLoss:
           tau_R = Losses_tau.get_ydata()
           dataToSave = np.append(dataToSave, [tau_R], axis=0)
        print 'Saving to file..'
        np.savetxt('RESULTS_DAT/%.20f.dat'%par.tt, dataToSave.T, header='%.7e %.7e %.7e'%(par.mu,par.g,par.R))
           
         
   plt.savefig('RESULTS/%.20f.png'%par.tt, bbox_inches='tight')
   plotCounter +=1

