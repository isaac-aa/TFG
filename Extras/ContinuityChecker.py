# The objective of this script is to check if for a given
# simulation, the mass is converved. It reads all the .dat in the
# given folder, and integrates numerically to obtain the total 'mass'

import glob
import numpy as np
import matplotlib.pyplot as plt

SimulationFolder = '../RESULTS/AnalyticalLossEq_HydrostaticP_SymV_NoDamping'



files = glob.glob(SimulationFolder+'/RESULTS_DAT/*.dat')

times = np.zeros(len(files))
integratedMass = np.zeros(len(files))

for i in range(len(times)):
    times[i] = float(files[i].split('/')[-1][:-4])
    z, rho = np.loadtxt(files[i], skiprows=0, usecols=(0,1), unpack=True)

    integratedMass[i] = np.trapz(rho[1:-1])

