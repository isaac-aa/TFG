import numpy as np
import matplotlib.pyplot as plt

R = 8.3144598e+07
mu = 1.113202
g = 27360.00
molarMass = 1.6605e-24

z_ref = 2.7555075e+07
z0 = 4.89736713*z_ref
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

def EulerPerturb(T0, dT):
   z_num = np.linspace(z[0], z[-1], 5000.)
   dz_num = z_num[1]-z_num[0]
   
   T_num = np.zeros(z_num.shape)
   T_num[0] = T0
   T_num[1] = T0 + dT*dz_num

   logp = np.zeros(z.shape)
   logp[0] = np.log(pA)
   logp[1] = logp[0] + dz*B/T[1]


   i = 2
   while i<len(T):
      rho = np.exp(logp[i-1])*mu / (R*T_num[i-1]) 
      numericalDensity = rho/(mu*molarMass)
      
      logLamda = np.interp(np.log10(T_num[i-1]), logT_table, logLamda_table)
      L_r = numericalDensity*numericalDensity*10**logLamda
      print dz_num*dz_num * 7.*L_r/(2.*A), T[i-1]**(7./2.)
     
      T7_2 = dz_num*dz_num * 7.*L_r/(2.*A) - T_num[i-2]**(7./2.) + 2*T[i-1]**(7./2.)
      T_num[i] = T7_2**(2./7.)

      logp[i] = logp[i-1] + dz*B/T[i]    

      i+=1
   
   return z_num, T_num, np.exp(logp)
   












