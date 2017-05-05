import Variables as var
import Flux

def LaxFriedichs():
   lamda = var.dt/var.dz

   massFlux, momentumFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentum, var.energy, var.v, var.P)

   var.rho[1:-1] = 0.5*(var.rho[2:]+var.rho[:-2]) - 0.5*lamda*(massFlux[2:] - massFlux[:-2]) 
   var.momentum[1:-1] = 0.5*(var.momentum[2:]+var.momentum[:-2]) - 0.5*lamda*(momentumFlux[2:] - momentumFlux[:-2]) + var.dt*var.momentumSource[1:-1]
   var.energy[1:-1] = 0.5*(var.energy[2:]+var.energy[:-2]) - 0.5*lamda*(energyFlux[2:] - energyFlux[:-2]) + var.dt*var.energySource[1:-1]







def FirstGen():
   lamda = var.dt/var.dz

   rhoHalf = var.rho
   momentumHalf = var.momentum
   energyHalf = var.energy


   # --------------------- LaxFriedichs Half Step
   # Compute change of variables

   massFlux, momentumFlux, energyFlux = Flux.ComputeFlux(var.rho, var.momentum, var.energy, var.v, var.P)

   rhoHalf[:-1] = 0.5*(var.rho[1:]+var.rho[:-1]) - 0.5*lamda*(massFlux[1:] - massFlux[:-1]) 
   momentumHalf[:-1] = 0.5*(var.momentum[1:]+var.momentum[:-1]) - 0.5*lamda*(momentumFlux[1:] - momentumFlux[:-1]) + 0.5*var.dt*var.momentumSource[:-1]
   energyHalf[:-1] = 0.5*(var.energy[1:]+var.energy[:-1]) - 0.5*lamda*(energyFlux[1:] - energyFlux[:-1]) + 0.5*var.dt*var.energySource[:-1]


   # Compute change of variables
   vHalf = momentumHalf/rhoHalf
   eHalf = energyHalf/rhoHalf - 0.5*vHalf*vHalf              #rho*e = E - 0.5*rho*v*v
   PHalf = rhoHalf*(var.gamma-1.)*eHalf


   massFluxHalf, momentumFluxHalf, energyFluxHalf = Flux.ComputeFlux(rhoHalf, momentumHalf, energyHalf, vHalf, PHalf)


   var.rho[1:-1] = var.rho[1:-1] - lamda*( massFluxHalf[1:-1] - massFluxHalf[:-2] )
   var.momentum[1:-1] = var.momentum[1:-1] - lamda*( momentumFluxHalf[1:-1] - momentumFluxHalf[:-2] ) + var.dt*var.momentumSource[1:-1]
   var.energy[1:-1] = var.energy[1:-1] - lamda*( energyFluxHalf[1:-1] - energyFluxHalf[:-2] ) + var.dt*var.energySource[1:-1]

   
