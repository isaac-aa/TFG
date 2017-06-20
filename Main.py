#--------------------------------------------
#		Main.py
# This module is the core of the simulation.
# Firstly it imports all the neccesary modules,
# creates the mesh, set the initial state and
# starts the main loop
#
#--------------------------------------------

import numpy as np
import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
comm_size = comm.Get_size()


CPUTime0 = time.clock()

import Parameters as par
import Grid

# ------------------ MESH CREATION -------------
Grid.Uniform1DGrid(par.N, par.z0, par.zf)
#Grid.ReadGridFromFile('Extras/ThermalLossesEq_IC.dat') #('hydrostatic_equilibrium_2.dat')

import Variables as var
import Settings as sets
import BoundaryConditions
import SourceTerm
import TimeStep
import ChangeOfVar
import Save
import SaveMPI

# ------------------ INITIAL CONDITION ---------

sets.InitialCondition(sets.argsIC)
ChangeOfVar.ConvertToPrim()
Save.Plot()
# ------------------ MAIN LOOP -----------------

if rank==0:
  CPUTimeConf = time.clock()
  print '################ \n'
  print 'Configuration took %.4f s'%(CPUTimeConf-CPUTime0)
  print '################ \t SIMULATION START \n'

comm.Barrier()

WallTime0 = time.time()
computeDTTime = 0
SchemeTime = 0
SourceTime = 0
BCTime = 0
SaveTime = 0

lastTime = 0
nowTime = 0

if comm_size > 1:
 if rank==0:
   sets.BoundaryConditionR = BoundaryConditions.InnerBoundary
   sets.argsR = np.zeros((4))
 elif rank==(comm_size-1):
   sets.BoundaryConditionL = BoundaryConditions.InnerBoundary
   sets.argsL = np.zeros((4))
 else:
   sets.BoundaryConditionR = BoundaryConditions.InnerBoundary
   sets.BoundaryConditionL = BoundaryConditions.InnerBoundary
   sets.argsL = np.zeros((4))
   sets.argsR = np.zeros((4))


while (par.it<=par.max_it and par.tt<=par.tf):
   #lastTime = time.clock()
   TimeStep.ComputeDT()  

   par.tt += par.dt  
   par.it += 1

   var.lastrho[:] = var.rho[:]
   var.lastmomentum[:] = var.momentum[:]
   var.lastenergy[:] = var.energy[:]
   #nowTime = time.clock()
   #computeDTTime += nowTime-lastTime

   # Time step
   #lastTime = nowTime
   sets.Scheme() 
   #nowTime = time.clock()
   #SchemeTime += nowTime-lastTime

   # Source computation
   #lastTime = nowTime
   if par.IsComputingSource:
      SourceTerm.ComputeSource()
   #nowTime = time.clock()
   #SourceTime += nowTime-lastTime

   # Update inner boundaries
   comm.Barrier()


   if comm_size > 1:
      buf = np.empty(4)
      # Now we are sending the information TO THE RIGHT (ie. recieving argsL) 
      if rank%2==0 and rank!=comm_size-1:
         comm.Send(np.array([0, var.rho[-2], var.momentum[-2], var.energy[-2]]), dest=rank+1, tag=1)
         #print 'Thread ', rank, ' sending [-2] to ', rank+1
      elif rank!=0:
         comm.Recv(buf, source=rank-1, tag=1)
         #print 'Thread ', rank, ' recieving ', buf[0], ' from ', rank-1
         sets.argsL[:] = buf[:]
      
      # Idem for odd subdomains
      if rank%2!=0 and rank!=comm_size-1:
         comm.Send(np.array([0, var.rho[-2], var.momentum[-2], var.energy[-2]]), dest=rank+1, tag=2)
         #print 'Thread ', rank, ' sending [-2] to ', rank+1
      elif rank!=0 and rank!=comm_size-1:
         comm.Recv(buf, source=rank-1, tag=2)
         sets.argsL[:] = buf[:]     
         #print 'Thread ', rank, ' recieving ', buf[0], ' from ', rank-1


      # Now we are sending the information TO THE LEFT (ie. recieving argsR)
      if rank%2==0 and rank!=0:
         comm.Send(np.array([1, var.rho[1], var.momentum[1], var.energy[1]]), dest=rank-1, tag=3)
         #print 'Thread ', rank, ' sending [1] to ', rank-1
      elif rank!=comm_size-1 and rank!=0:
         comm.Recv(buf, source=rank+1, tag=3)
         sets.argsR[:] = buf[:]
         #print 'Thread ', rank, ' recieving ', buf[0], ' from ', rank+1
      
      # Idem for odd subdomains
      if rank%2!=0 and rank!=0:
         comm.Send(np.array([1, var.rho[1], var.momentum[1], var.energy[1]]), dest=rank-1, tag=4)
         #print 'Thread ', rank, ' sending [1] to ', rank-1
      elif rank!=comm_size-1:
         comm.Recv(buf, source=rank+1, tag=4)
         #print 'Thread ', rank, ' recieving ', buf[0], ' from ', rank+1
         sets.argsR[:] = buf[:]
      

   comm.Barrier()
   # Boundary conditions     
   #lastTime = nowTime
   sets.BoundaryConditionL(sets.argsL)
   sets.BoundaryConditionR(sets.argsR)
   #nowTime = time.clock()
   #BCTime = nowTime-lastTime

   # Compute change of variables
   ChangeOfVar.ConvertToPrim()

   if comm_size == 0:
    if par.it%par.save_rate == 0.:
     
     #lastTime = time.clock()
     Save.Plot()
     #nowTime = time.clock()
     #SaveTime += nowTime-lastTime

     print '\n   ## IT: %d \t CPU-Time: %.2f s \t Wall-time: %.2f'%(par.it, nowTime-CPUTimeConf, time.time()-WallTime0)
     print 'DT: %.3e \t t: %.3f s \t CFL: %.3e'%(par.dt, par.tt, par.cfl)

     #nowTime /= 100.
     #print ' ComputeDT \t %.2f \n Scheme \t %.2f \n Source \t %.2f \n BC \t\t %.2f \n Save \t\t %.2f'%(computeDTTime/nowTime, SchemeTime/nowTime, SourceTime/nowTime, BCTime/nowTime, SaveTime/nowTime)


     #maxIndex = np.argmax(np.abs(var.v))
     #print 'Max v: %.2f cm/s @ z[%d]: %.3e \n'%(var.v[maxIndex], maxIndex, Grid.z[maxIndex])
   else:
     if par.it%par.save_rate == 0.:
       SaveMPI.CollectiveSave()
      
       if rank==0:
         nowTime = time.clock()
         print '\n   ## IT: %d \t CPU-Time: %.2f s \t Wall-time: %.2f'%(par.it, nowTime-CPUTimeConf, time.time()-WallTime0)
         print 'DT: %.3e \t t: %.3f s \t CFL: %.3e'%(par.dt, par.tt, par.cfl)



if rank==0:

   print '################ \t SIMULATION ENDED \n'
   CPUTimeEnd = time.clock()
   print 'Simulation took %.4f s'%(CPUTimeEnd-CPUTime0)
   print '################ \n'











