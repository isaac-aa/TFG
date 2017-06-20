#--------------------------------------------
#		SaveMPI.py
# This module saves the current state as a file,
# which is handled using MPI and collective I/O
# This way, all proccesses can write to the same
# file without obstructions
#
#--------------------------------------------


import shutil
import os
import numpy as np
import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

amode = MPI.MODE_WRONLY|MPI.MODE_CREATE

import Grid
import Parameters as par
import Variables as var
import Settings as sets
import Analytic
import Characteristics

buf_size = Grid.N_sub-2  # The BC are not stored
buffer_rho = np.empty(buf_size, dtype='float') 
buffer_mom = np.empty(buf_size, dtype='float')
buffer_ene = np.empty(buf_size, dtype='float')

def CollectiveSave():
    """
    Saves to a single file for each dt all the data stored in the proccesses
    """
    filename = par.FolderName + '/RESULTS_DAT/%.20f.dat'%par.tt 
    fh = MPI.File.Open(comm, filename, amode)

    buffer_rho = var.rho[1:-1]
    buffer_mom = var.momentum[1:-1]
    buffer_ene = var.energy[1:-1]
    
    offset_rho = rank*buffer_rho.nbytes
    offset_mom = rank*buffer_mom.nbytes + buffer_mom.nbytes*size
    offset_ene = rank*buffer_ene.nbytes + 2*buffer_mom.nbytes*size
    
    fh.Write_at_all(offset_rho, buffer_rho)
    fh.Write_at_all(offset_mom, buffer_mom)
    fh.Write_at_all(offset_ene, buffer_ene)
       
    fh.Close() 
