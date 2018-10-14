#!/usr/bin/env python
from __future__ import print_function, division
import sys
import numpy as np
import MDAnalysis as mda
import matplotlib
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms
import time
from shutil import copyfile
import glob, os
from MDAnalysis import Writer
import mpi4py
from mpi4py import MPI

#---------------------------------------
MPI.Init

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
if rank == 0:
    print(mda.__version__)

#------------------------------------------
j = sys.argv[1]

def block_rmsd(index, topology, trajectory, xref0, start=None, stop=None, step=None):
    start00 = time.time()
    clone = mda.Universe(topology, trajectory)
    g = clone.atoms[index]
    
    print("block_rmsd", start, stop, step)
    
    bsize = int(stop-start)
    results = np.zeros([bsize,2], dtype=float)
    t_comp = np.zeros(bsize, dtype=float)
    t_IO = np.zeros(bsize, dtype=float)

    start1 = time.time()
    start0 = start1
    for iframe, ts in enumerate(clone.trajectory[start:stop:step]):
        start2 = time.time()
        results[iframe, :] = ts.time, rms.rmsd(g.positions, xref0, center=True, superposition=True)
        t_comp[iframe] = time.time()-start2
        t_IO[iframe] = start2-start1
        start1 = time.time()
    
    t_end_loop = time.time()-start1
    t_all_frame = time.time()-start0
    t_IO_final = np.sum(t_IO)
    t_comp_final = np.sum(t_comp)  
    t_init = start0-start00
    return results, t_comp_final, t_IO_final, t_all_frame, t_end_loop, t_init 
#-----------------------------------------------------------------------
# Check the files in the directory
PSF = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'files/adk4AKE.psf')))
filenames = os.listdir(os.path.join(os.getcwd(),'files'))
longXTC = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'newtraj.xtc')))
longXTC1 = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'files/newtraj{}.xtc'.format(j))))

if rank == 0:
   copyfile(longXTC, longXTC1)
   u = mda.Universe(PSF, longXTC1)
   print(len(u.trajectory))

MPI.COMM_WORLD.Barrier()

start1 = time.time()
#----------------------------------------------------------------------
u = mda.Universe(PSF, longXTC1)
mobile = u.select_atoms("(resid 1:29 or resid 60:121 or resid 160:214) and name CA")[1:147]
index = mobile.indices
topology, trajectory = mobile.universe.filename, mobile.universe.trajectory.filename
ref0 = mobile
xref0 = ref0.positions-ref0.center_of_mass()

# Create each segment for each process
start2 = time.time()
   
frames_seg = np.zeros([size,2], dtype=int)
bsize = int(np.ceil(mobile.universe.trajectory.n_frames / float(size)))
for iblock in range(size):
    frames_seg[iblock, :] = iblock*bsize, (iblock+1)*bsize

d = dict([key, frames_seg[key]] for key in range(size))  

start, stop = d[rank][0], d[rank][1] 

# Block-RMSD in Parallel
start3 = time.time()
out = block_rmsd(index, topology, trajectory, xref0, start=start, stop=stop, step=1) 

# Communication
start4 = time.time()
if rank == 0:
   data1 = np.zeros([size*bsize,2], dtype=float)
else:
   data1 = None

comm.Gather(out[0], data1, root=0)

start5 = time.time()
if rank == 0:
   data = np.zeros([size,5], dtype=float)
else:
   data = None

comm.Gather(np.array(out[1:], dtype=float), data, root=0)

start6 = time.time()

# Cost Calculation
init_time = start2-start1
comm_time1 = start3-start2
comm_time2 = start5-start4
comm_time3 = start6-start5
comp_time = start4-start3
tot_time = comp_time+comm_time2

tot_time = comm.gather(tot_time, root=0)
init_time = comm.gather(init_time, root=0)
comm_time1 = comm.gather(comm_time1, root=0)
comm_time2 = comm.gather(comm_time2, root=0)
comm_time3 = comm.gather(comm_time3, root=0)
comp_time = comm.gather(comp_time, root=0)

if rank == 0:
   tot_time = tot_time
   init_time = init_time
   comm_time1 = comm_time1
   comm_time2 = comm_time2
   comm_time3 = comm_time3
   comp_time = comp_time
else:
   tot_time = None
   init_time = None
   comm_time1 = None
   comm_time2 = None
   comm_time3 = None
   comp_time = None


# Storing the data
if rank == 0:
   os.remove('files/newtraj{}.xtc'.format(j))   
   os.remove('files/.newtraj{}.xtc_offsets.npz'.format(j))
   print(tot_time, comm_time1, comm_time2)
   with open('data1.txt', mode='a') as file1:
        for i in range(size):
            file1.write("{} {} {} {} {} {} {} {}\n".format(size, j, i, data[i][0], data[i][1], data[i][2], data[i][3], data[i][4]))
   with open('data2.txt', mode='a') as file2:
        file2.write("{} {} {}\n".format(size, j, tot_time))
        file2.write("{} {} {}\n".format(size, j, init_time))
        file2.write("{} {} {}\n".format(size, j, comm_time1)) 
        file2.write("{} {} {}\n".format(size, j, comm_time2))
        file2.write("{} {} {}\n".format(size, j, comm_time3))
        file2.write("{} {} {}\n".format(size, j, comp_time))	

MPI.Finalize

