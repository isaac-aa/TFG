import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

R = 8.3144598e+07
mu = 1.113202
g = 27360.00
molarMass = 1.6605e-24

z_ref = 2.7555075e+07
z0 = 144352725.54647857 #4.89736713*z_ref
zf = 33.33594543*z_ref
z = np.linspace(z0, zf, 5000)

dz = z[1]- z[0]

TA = 1e5
TB = 1e6
rhoA = 0.06589083*5.7181381e-13
pA = rhoA*R*TA/mu 

A = 9e-12 * 1e5
B = -mu*g/R

D = (zf -z0*(TB/TA)**(7./2.))/( (TB/TA)**(7./2.) -1. )
alpha = TA**(7./2.)/(z0+D)
logp0 = np.log(pA) - B*alpha**(-2./7.)*(z0+D)**(5./7.)*7./5.



T = (alpha*(z+D))**(2./7.)
logp = logp0 + B*alpha**(-2./7.)*(z+D)**(5./7.)*7./5.
p = np.exp(logp)

rho = p*mu/(R*T)

def Euler():
   T = np.zeros(z.shape)
   T[0] = TA
   
   logp = np.zeros(z.shape)
   logp[0] = np.log(pA)

   i=1
   while i<len(T):
      T[i] = T[i-1] +  dz*T[i-1]**(-5./2.)*alpha*2./7.
      logp[i] = logp[i-1] + dz*B/T[i]
      
      i+=1 

   return T, np.exp(logp)

dT = (T[1]-T[0])/dz

logT_table, Lamda_table = np.loadtxt('../dere_etal_table.dat', usecols=(0,1), unpack=True)
logLamda_table = np.log10(Lamda_table)

#dT = .279 buen valor
def EulerPerturb(rhoA, T0, dT):
   z_num = np.linspace(z[0], z[-1], 5000)
   dz_num = z_num[1]-z_num[0]
   
   T_num = np.zeros(z_num.shape)
   T_num[0] = T0
   T_num[1] = T0 + dT*dz_num
 
   pA = rhoA*R*TA/mu

   logp = np.zeros(z_num.shape)
   logp[0] = np.log(pA)
   logp[1] = logp[0] + dz*B/T_num[1]


   i = 2
   L_r = np.zeros(z_num.shape)
   while i<len(z_num):
      rho = np.exp(logp[i-1])*mu / (R*T_num[i-1]) 
      numericalDensity = rho/(mu*molarMass)
      
      logLamda = np.interp(np.log10(T_num[i-1]), logT_table, logLamda_table)
      
      L_r[i] = numericalDensity*numericalDensity*10**logLamda
      # -dz_num*dz_num*7.*L_r[i]/(2.*A)
      # Parece ir con '+'. Sin embargo, las soluciones o no alcanzan 1e6K o lo hacen
      # de manera demasiado similar a el caso Lr=0
      T7_2 = - T_num[i-2]**(7./2.) + 2.*T_num[i-1]**(7./2.)  -  dz_num*dz_num*7.*L_r[i]/(2.*A)
      T_num[i] = T7_2**(2./7.) #+ ( dz_num*dz_num*7.*L_r[i]/(2.*A) )

      logp[i] = logp[i-1] + dz_num*B/T_num[i]    

      i+=1
   
   kappa = A*T_num**(5./2.)
   dT = (T_num[2:]-T_num[:-2])/dz_num
   F = -kappa[1:-1]*dT
   dF = (F[2:]-F[:-2])/dz_num
   return z_num, T_num, np.exp(logp), L_r, F, dF
   


dT2 = (T[-2]-T[-1])/dz  #Buen valor 0.00025 o 0.00013
def EulerPerturbInvert(T0, dT):
   z_num = np.linspace(z[0], z[-1], 5000.)
   dz_num = np.abs(z_num[1]-z_num[0])
   
   T_num = np.zeros(z_num.shape)
   T_num[-1] = T0
   T_num[-2] = T0 + dT*dz_num

   logp = np.zeros(z_num.shape)
   logp[-1] = np.log(0.20231163238863872) #np.log(p[-1])   
   logp[-2] = logp[-1] - dz_num*B/T_num[-2]


   i = -3
   while abs(i)<=len(z_num):
      rho = np.exp(logp[i+1])*mu / (R*T_num[i+1]) 
      numericalDensity = rho/(mu*molarMass)
      
      logLamda = np.interp(np.log10(T_num[i+1]), logT_table, logLamda_table)
      
      L_r = numericalDensity*numericalDensity*10**logLamda
      #print L_r
      #print (dz_num*dz_num * 7.*L_r/(2.*A) )**(2./7.), T[i-1]
      #L_r = 0.

      T7_2 = -dz_num*dz_num*7.*L_r/(2.*A) - T_num[i+2]**(7./2.) + 2*T_num[i+1]**(7./2.)
      T_num[i] = T7_2**(2./7.) 
      #print T_num[i+2], T_num[i+1], T_num[i]
      logp[i] = logp[i+1] - dz_num*B/T_num[i]    

      i-=1
   
   return z_num, T_num, np.exp(logp)


#z_num1, T_num1, p_num1 = EulerPerturb(T[0], .278)
#rho1 = p_num1*mu/(R*T_num1)

#np.savetxt('ThermalLossesEq_ana.dat', np.array([z_num1, rho1, p_num1, T_num1]).T )


z_in = z-dz/2.
z_in = np.append(z_in, z[-1]+dz/2.)

p_in = np.interp(z_in, z, p)

p_in[0] = 2*p[0] - p_in[1]
p_in[-1] = 2*p[-1] - p_in[-2]

T_in = np.interp(z_in, z, T)

T_in[0] = 2*T[0] - T_in[1]
T_in[-1] = 2*T[-1] - T_in[-2]

rho_in = np.interp(z_in, z, rho)

rho_in[0] = 2*rho[0] - rho_in[1]
rho_in[-1] = 2*rho[-1] - rho_in[-2]

#np.savetxt('ThermalEq_IC.dat', np.array([z_in, p_in, rho_in]).T )
#np.savetxt('ThermalEq_IC_2.dat', np.array([z, p, rho]).T )












