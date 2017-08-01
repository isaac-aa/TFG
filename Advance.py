#--------------------------------------------
#		Advance.py
# This module defines the different numerical
# schemes that can be used in this code. The
# scheme should be set at Settings.py
#
#--------------------------------------------



import Parameters as par
import Grid
import Variables as var
import numpy as np
import Flux

def LaxFriedrichs():
   lamda = par.dt/Grid.dz

   massFlux, momentumFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentum, var.energy, var.v, var.P)

   var.rho[1:-1] = 0.5*(var.rho[2:]+var.rho[:-2]) - 0.5*lamda*(massFlux[2:] - massFlux[:-2]) 
   var.momentum[1:-1] = 0.5*(var.momentum[2:]+var.momentum[:-2]) - 0.5*lamda*(momentumFlux[2:] - momentumFlux[:-2]) 
   var.energy[1:-1] = 0.5*(var.energy[2:]+var.energy[:-2]) - 0.5*lamda*(energyFlux[2:] - energyFlux[:-2]) 




def LaxFriedrichs2D():
   lamdaZ = par.dt/Grid.dz  # Minus because for np.roll(-, 1, axis=1) means the i-1 value
   lamdaY = par.dt/Grid.dy

   massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY = Flux.ComputeFlux2D(var.rho, var.momentumZ, var.momentumY, var.energy, var.vZ, var.vY, var.P)

   var.rho[1:-1, 1:-1] = 0.25*(var.rho[2:, 1:-1] + var.rho[:-2, 1:-1] + var.rho[1:-1, 2:] + var.rho[1:-1, :-2]) - 0.5*lamdaZ*(massFluxZ[1:-1, 2:] - massFluxZ[1:-1, :-2]) - 0.5*lamdaY*(massFluxY[2:, 1:-1] - massFluxY[:-2, 1:-1])

   var.momentumZ[1:-1, 1:-1] = 0.25*(var.momentumZ[2:, 1:-1] + var.momentumZ[:-2, 1:-1] + var.momentumZ[1:-1, 2:] + var.momentumZ[1:-1, :-2]) - 0.5*lamdaZ*(momentumZFluxZ[1:-1, 2:] - momentumZFluxZ[1:-1, :-2]) - 0.5*lamdaY*(momentumZFluxY[2:, 1:-1] - momentumZFluxY[:-2, 1:-1])

   var.momentumY[1:-1, 1:-1] = 0.25*(var.momentumY[2:, 1:-1] + var.momentumY[:-2, 1:-1] + var.momentumY[1:-1, 2:] + var.momentumY[1:-1, :-2]) - 0.5*lamdaY*(momentumYFluxZ[1:-1, 2:] - momentumYFluxZ[1:-1, :-2]) - 0.5*lamdaY*(momentumYFluxY[2:, 1:-1] - momentumYFluxY[:-2, 1:-1]) 

   var.energy[1:-1, 1:-1] = 0.25*(var.energy[2:, 1:-1] + var.energy[:-2, 1:-1] + var.energy[1:-1, 2:] + var.energy[1:-1, :-2]) - 0.5*lamdaZ*(energyFluxZ[1:-1, 2:] - energyFluxZ[1:-1, :-2]) - 0.5*lamdaY*(energyFluxY[2:, 1:-1] - energyFluxY[:-2, 1:-1])


def FirstGen():
   lamda = par.dt/Grid.dz

   rhoHalf = np.ones(Grid.z.shape)
   momentumHalf = np.ones(Grid.z.shape)
   energyHalf = np.ones(Grid.z.shape)


   # --------------------- LaxFriedichs Half Step
   # Compute change of variables

   massFlux, momentumFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentum, var.energy, var.v, var.P)

   rhoHalf[:-1] = 0.5*(var.rho[1:]+var.rho[:-1]) - 0.5*lamda*(massFlux[1:] - massFlux[:-1]) 
   momentumHalf[:-1] = 0.5*(var.momentum[1:]+var.momentum[:-1]) - 0.5*lamda*(momentumFlux[1:] - momentumFlux[:-1]) 
   energyHalf[:-1] = 0.5*(var.energy[1:]+var.energy[:-1]) - 0.5*lamda*(energyFlux[1:] - energyFlux[:-1]) 


   # Compute change of variables
   vHalf = momentumHalf/rhoHalf
   eHalf = energyHalf/rhoHalf - 0.5*vHalf*vHalf              #rho*e = E - 0.5*rho*v*v
   PHalf = rhoHalf*(par.gamma-1.)*eHalf


   massFluxHalf, momentumFluxHalf, energyFluxHalf = Flux.ComputeFlux(rhoHalf, momentumHalf, energyHalf, vHalf, PHalf)


   var.rho[1:-1] = var.rho[1:-1] - lamda*( massFluxHalf[1:-1] - massFluxHalf[:-2] )
   var.momentum[1:-1] = var.momentum[1:-1] - lamda*( momentumFluxHalf[1:-1] - momentumFluxHalf[:-2] ) 
   var.energy[1:-1] = var.energy[1:-1] - lamda*( energyFluxHalf[1:-1] - energyFluxHalf[:-2] )




def LaxWendroffRitchmyer2D():
   lamdaZ = par.dt/Grid.dz
   lamdaY = par.dt/Grid.dy

   rhoHalf = np.zeros(Grid.z.shape)
   momentumZHalf = np.zeros(Grid.z.shape)
   momentumYHalf = np.zeros(Grid.z.shape)
   energyHalf = np.zeros(Grid.z.shape)


   # --------------------- LaxFriedichs Half Step
   massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY = Flux.ComputeFlux2D(var.rho, var.momentumZ, var.momentumY, var.energy, var.vZ, var.vY, var.P)
 
   # CHECK REPEATED INDEX IN THE MEAN COMPUTATION
   rhoHalf[:-1, :-1] = 0.25*(var.rho[1:,:-1]+var.rho[:-1,:-1] + var.rho[:-1,1:]+ var.rho[:-1,:-1]) - 0.5*lamdaZ*(massFluxZ[:-1,1:] - massFluxZ[:-1,:-1]) - 0.5*lamdaY*(massFluxY[1:, :-1] - massFluxY[:-1, :-1]) 
  
   momentumZHalf[:-1, :-1] = 0.25*(var.momentumZ[1:,:-1]+var.momentumZ[:-1,:-1] + var.momentumZ[:-1,1:] + var.momentumZ[:-1,:-1]) - 0.5*lamdaZ*(momentumZFluxZ[:-1, 1:] - momentumZFluxZ[:-1,:-1]) - 0.5*lamdaY*(momentumZFluxY[1:,:-1]-momentumZFluxY[:-1,:-1]) 
   
   momentumYHalf[:-1, :-1] = 0.25*(var.momentumY[1:,:-1]+var.momentumY[:-1,:-1] + var.momentumY[:-1,1:] + var.momentumY[:-1,:-1]) - 0.5*lamdaZ*(momentumYFluxZ[:-1, 1:] - momentumYFluxZ[:-1,:-1]) - 0.5*lamdaY*(momentumYFluxY[1:,:-1]-momentumYFluxY[:-1,:-1]) 

   energyHalf[:-1,:-1] = 0.25*(var.energy[1:,:-1]+var.energy[:-1,:-1] + var.energy[:-1,1:]+var.energy[:-1,:-1]) - 0.5*lamdaZ*(energyFluxZ[:-1,1:] - energyFluxZ[:-1,:-1]) - 0.5*lamdaY*(energyFluxY[1:,:-1]-energyFluxY[:-1,:-1]) 


   # Compute change of variables
   vZHalf = momentumZHalf/rhoHalf
   vYHalf = momentumYHalf/rhoHalf
   eHalf = energyHalf/rhoHalf - 0.5*(vZHalf*vZHalf + vYHalf*vYHalf)              #rho*e = E - 0.5*rho*v*v
   PHalf = rhoHalf*(par.gamma-1.)*eHalf

   massFluxZ, massFluxY, momentumZFluxZ, momentumZFluxY, momentumYFluxZ, momentumYFluxY, energyFluxZ, energyFluxY = Flux.ComputeFlux2D(rhoHalf, momentumZHalf, momentumYHalf, energyHalf, vZHalf, vYHalf, PHalf)


   var.rho[1:-1, 1:-1] = var.rho[1:-1, 1:-1] - lamdaZ*( massFluxZ[1:-1, 1:-1] - massFluxZ[1:-1, :-2] ) - lamdaY*( massFluxY[1:-1, 1:-1] - massFluxY[:-2, 1:-1] )
   var.momentumZ[1:-1, 1:-1] = var.momentumZ[1:-1, 1:-1] - lamdaZ*( momentumZFluxZ[1:-1, 1:-1] - momentumZFluxZ[1:-1,:-2] ) - lamdaY*( momentumZFluxY[1:-1,1:-1]- momentumZFluxY[:-2,1:-1]) 
   var.momentumY[1:-1, 1:-1] = var.momentumY[1:-1, 1:-1] - lamdaZ*( momentumYFluxZ[1:-1, 1:-1] - momentumYFluxZ[1:-1,:-2] ) - lamdaY*( momentumYFluxY[1:-1,1:-1]- momentumYFluxY[:-2,1:-1]) 
   var.energy[1:-1, 1:-1] = var.energy[1:-1, 1:-1] - lamdaZ*( energyFluxZ[1:-1, 1:-1] - energyFluxZ[1:-1,:-2] ) - lamdaY*( energyFluxY[1:-1,1:-1] - energyFluxY[:-2,1:-1] )

   


   
