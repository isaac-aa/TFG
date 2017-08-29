#--------------------------------------------
#		Main.py
# This module is the core of the simulation.
# Firstly it imports all the neccesary modules,
# creates the mesh, set the initial state and
# starts the main loop
#
#--------------------------------------------

import time

CPUTime0 = time.clock()

import Parameters as par
import Grid

# ------------------ MESH CREATION -------------
Grid.Uniform2DGrid(par.Nz, par.Ny, par.z0, par.zf, par.y0, par.yf)
#Grid.Uniform1DGrid(par.Nz, par.z0, par.zf)
#Grid.ReadGridFromFile('Extras/ThermalEq_IC_2.dat') #('hydrostatic_equilibrium_2.dat')

import Settings as sets
import SourceTerm
import TimeStep
import ChangeOfVar
import Save
import Variables as var

# ------------------ INITIAL CONDITION ---------

sets.InitialCondition(sets.argsIC)
ChangeOfVar.ConvertToPrim()

# ------------------ MAIN LOOP -----------------

Save.Plot()




CPUTimeConf = time.clock()
print '################ \n'
print 'Configuration took %.4f s'%(CPUTimeConf-CPUTime0)


# ------------------ MAIN LOOP -----------------

print '################ \t SIMULATION START \n'

WallTime0 = time.time()
computeDTTime = 0
SchemeTime = 0
SourceTime = 0
BCTime = 0
SaveTime = 0

lastTime = 0
nowTime = 0
while (par.it<=par.max_it and par.tt<=par.tf):
   lastTime = time.clock()
   TimeStep.ComputeDT()  
   var.updateLastVars()
   par.tt += par.dt  
   par.it += 1
   nowTime = time.clock()
   computeDTTime += nowTime-lastTime
   
   # Time step
   lastTime = nowTime
   sets.Scheme.compute()
   nowTime = time.clock()
   SchemeTime += nowTime-lastTime
   ChangeOfVar.ConvertToPrim()
   
   
   #SourceTerm.computeImplicitConduction()   
   #"""
   #Source computation
   lastTime = nowTime
   if par.IsComputingSource:
      SourceTerm.ComputeSource()
   nowTime = time.clock()
   SourceTime += nowTime-lastTime
   #"""

   # Boundary conditions     
   lastTime = nowTime
   if sets.BoundaryConditionL!=None: sets.BoundaryConditionL.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
   if sets.BoundaryConditionR!=None: sets.BoundaryConditionR.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
   if par.dim==2:
      if sets.BoundaryConditionT != None: sets.BoundaryConditionT.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
      if sets.BoundaryConditionB != None: sets.BoundaryConditionB.computeBC(var.rho, (var.momentumZ, var.momentumY), var.energy)
   
   
   nowTime = time.clock()
   BCTime = nowTime-lastTime
   
   # Compute change of variables
   ChangeOfVar.ConvertToPrimBoundaries()
   
   if sets.Tracers != None:
      for Tracer in sets.Tracers:
         Tracer.TimeStep()

   if par.it%par.save_rate == 0.:
     

      lastTime = time.clock()
      Save.Plot()
      nowTime = time.clock()
      SaveTime += nowTime-lastTime


      print '\n   ## IT: %d \t CPU-Time: %.2f s \t Wall-time: %.2f'%(par.it, nowTime-CPUTimeConf, time.time()-WallTime0)
      print 'DT: %.3e \t t: %.3f s \t CFL: %.3e'%(par.dt, par.tt, par.cfl)

      nowTime /= 100.
      print ' ComputeDT \t %.2f \n Scheme \t %.2f \n Source \t %.2f \n BC \t\t %.2f \n Save \t\t %.2f'%(computeDTTime/nowTime, SchemeTime/nowTime, SourceTime/nowTime, BCTime/nowTime, SaveTime/nowTime)









