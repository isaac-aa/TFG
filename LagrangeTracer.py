import Grid
import Variables as var
import Parameters as par
import numpy as np

class Tracer():
   z = 0.
   vZ = 0.
   vY = 0.
   
   def __init__(self, z, vz):
      self.z = z
      self.vZ = vz
   
   def TimeStep(self):
      self.z += par.dt*np.interp(self.z, Grid.z, var.vZ)