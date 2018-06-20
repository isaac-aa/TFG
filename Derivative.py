"""
Module to compute the derivative of the variables at the domain
It follows a object structure and must be set at Settings.py

@author: Isaac Alonso
"""
import Grid

class GenericDerivative(object):
   name = 'Generic Derivative'
   dim = 0

   def __str__(self):
      return 'Derivative scheme: ' + self.name + ' ; dim=', str(self.dim) 


   def setup(self):
      """
      Setup function to pre-allocate needed data
      """
      print 'Setup \t '+ str(self.__str__())


class CentralDer2D(GenericDerivative):
   """
   Second order centered derivative in two dimensions
   """
   name = 'Second order centered'
   dim = 2
   def __call__(self,variable, axis):
      if axis==0:  # z-axis
         return 0.5*(variable[1:-1,2:]-variable[1:-1,:-2])/Grid.dz
      elif axis==1:  # y-axis
         return 0.5*(variable[2:,1:-1]-variable[:-2,1:-1])/Grid.dy
      else:
         print '## ERROR ## \t Axis not defined for ', self.name
         exit()

