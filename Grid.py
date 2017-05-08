#--------------------------------------------
#		Grid.py
# This module creates the mesh for the simulation, 
# using as inputs the parameters given by the user
# at Parameters.py
#
#--------------------------------------------

import numpy as np

import Parameters as par

z = np.array([])  #Empty array that will be filled with the grid data
dz = 0.

def Uniform1DGrid(Ncell, z0, zf):
   global z
   global dz
   
   dz = (zf-z0)/Ncell
   z = np.arange(z0-dz/2., zf+3*dz/2., dz)  #If N=100, this will yield N+2 nodes (ghost nodes)




