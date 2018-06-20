#--------------------------------------------
#		Save.py
# This module saves the current state as a plot,
# or a ASCII file (TODO). It can add different
# plots, that can be set at Parameters.py (eg.
# an analytic solution)
#
#--------------------------------------------

import shutil
import os
import numpy as np
import time

import Grid
import Parameters as par
import Variables as var
import Settings as sets


# ------------------ PLOT ----------------------


print 'Loading Save..'

# Create results folders
if not os.path.exists(par.FolderName):
   os.makedirs(par.FolderName)
   os.makedirs(par.FolderName + '/RESULTS')
   os.makedirs(par.FolderName + '/RESULTS_DAT')
  

ThereIsSettings = os.path.exists(par.FolderName+'/Settings.py')
ThereIsParameters = os.path.exists(par.FolderName+'/Parameters.py')

if ThereIsSettings or ThereIsParameters:
   print '###### WARNING ######'
   print 'There is already a previous simulation stored at ' + par.FolderName
   des = raw_input('Do you want to overwrite it? ([y]/n)')
   fh = open(par.FolderName + '/log.txt', 'w')
   fh.write('#Starting simulation\n')
   fh.close()
   if des ==  'n':
      exit()

shutil.copy2('Settings.py', par.FolderName)
shutil.copy2('Parameters.py', par.FolderName)

ZeroTime = time.clock()
plotCounter = 0

def Plot():
   
   np.save(par.FolderName + '/RESULTS_DAT/%d.dat'%plotCounter, np.array([Grid.z, Grid.y, var.rho, var.momentumZ, var.momentumY, var.energy, var.momentumX,var.Bz, var.By, var.Bx]) ) 
      
   fh = open(par.FolderName + '/log.txt','a')
   fh.write(str(par.it) + ' ' + str(par.dt) + ' ' + str(par.tt) + ' ' + str(time.clock() - ZeroTime) + '\n')     
   fh.close()
   
   
   if sets.Tracers != None:
      fh = open(par.FolderName + '/RESULTS_DAT/%.20f.tcr'%par.tt,'w')
      for Trace in sets.Tracers:
         fh.write(str(Trace.z) + ' ' + str(Trace.vZ) + '\n')     
      fh.close()
      




   
   





