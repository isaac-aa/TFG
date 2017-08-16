#--------------------------------------------
#		Grid.py
# This module creates the mesh for the simulation, 
# using as inputs the parameters given by the user
# at Parameters.py
#
#--------------------------------------------

import numpy as np

import Parameters as par

z = np.array([])  # Empty array that will be filled with the grid data
y = np.array([])  # Same shape, but for Y-axis coordinates
dz = 0.
dy = 0.

def Uniform1DGrid(Ncell, z0, zf):
   global z
   global dz
   
   dz = (zf-z0)/Ncell
   z = np.arange(z0-dz/2., zf+3*dz/2., dz)  #If N=100, this will yield N+2 nodes (ghost nodes)

def Uniform2DGrid(NcellZ, NcellY, z0, zf, y0, yf):
   global z
   global y
   global dz
   global dy
   
   dz = (zf-z0)/NcellZ
   dy = (yf-y0)/NcellY

   z_ax = np.arange(z0-dz/2, zf+3*dz/2., dz)
   y_ax = np.arange(y0-dy/2, yf+3*dy/2., dy)

   z, y = np.meshgrid(z_ax, y_ax)



def ReadGridFromFile(FileName):  
   print "Creating mesh from " + FileName
   global z
   global dz
   
   f = open(FileName)
   par.N = int(f.readline().split()[0])
   z_ref = float(f.readline().split()[0])
   f.close()
   

   data = np.loadtxt(FileName, skiprows=2)
   data = data.T
   z = data[0]*z_ref
   dz = z[1]-z[0]
   print "Mesh created with Nz=%d and dz=%.3e"%(par.Nz, dz)
   
   
   
   
   
   

   
