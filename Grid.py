#--------------------------------------------
#		Grid.py
# This module creates the mesh for the simulation, 
# using as inputs the parameters given by the user
# at Parameters.py
#
#--------------------------------------------

import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()




import Parameters as par

z = np.array([])  #Empty array that will be filled with the grid data
dz = 0.

N_sub = 0 # Number of nodes per subdomain

def Uniform1DGrid(Ncell, z0, zf):
   global z
   global dz
   global N_sub

   subdomains = comm.Get_size()
  

   dz_i = (zf-z0)/subdomains  # Space for each node
   dz = (zf-z0)/Ncell         # Spatial resolution
   
   z = np.arange( (z0+rank*dz_i)-dz/2, (z0 + (rank+1)*dz_i)+3*dz/2, dz) 
   N_sub = len(z)

   # Old code for one node
   #dz = (zf-z0)/Ncell
   #z = np.arange(z0-dz/2., zf+3*dz/2., dz)  #If N=100, this will yield N+2 nodes (ghost nodes)

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
   print "Mesh created with N=%d and dz=%.3e"%(par.N, dz)
   
   
   


   
   
   

   
