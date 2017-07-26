import numpy as np
import shutil
import os

import Grid
import Parameters as par
import Variables as var
import Settings as sets



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
   if des ==  'n':
      exit()

shutil.copy2('Settings.py', par.FolderName)
shutil.copy2('Parameters.py', par.FolderName)

def Plot2D():
   """
   Saves the current state to a binary file.
   The information is stored in rows for each node, where the colums are the stored variables, namely: z, y, rho, momentumZ, momentumY, energy
   """

   np.save(par.FolderName + '/RESULTS_DAT/%.20f.dat'%par.tt, np.array([Grid.z, Grid.y, var.rho, var.momentumZ, var.momentumY, var.energy])) 
