#!/usr/bin/env python
from __future__ import print_function, division
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import time
from shutil import copyfile
import glob, os
from MDAnalysis import Writer
import mpi4py
from mpi4py import MPI
import h5py 

#---------------------------------------
MPI.Init

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
if rank == 0:
   print (mda.__version__)

#---------------------------------------
j = sys.argv[1]

def block_rmsd(dset, xref0, start=None, stop=None, step=None):
    start00 = time.time()
   
    print("block_rmsd", start, stop, step)

    bsize = int(stop-start)
    results = np.zeros([bsize,2], dtype=float)
    t_comp = np.zeros(bsize, dtype=float)
    t_IO = np.zeros(bsize, dtype=float)

    start0 = time.time()
    for iframe, ts in enumerate(range(start,stop,step)):
        start1 = time.time()
        pos = np.reshape(dset[ts,:], (-1,3))
        start2 = time.time()
        results[iframe, :] = ts, rms.rmsd(pos, xref0, center=True, superposition=True)
        t_comp[iframe] = time.time()-start2
        t_IO[iframe] = start2-start1
        start1 = time.time()

    start3 = time.time()
    t_end_loop = start3-start1
    t_all_frame = start3-start0
    t_IO_final = np.sum(t_IO)
    t_comp_final = np.sum(t_comp)
    t_init = start0-start00
    return results, t_comp_final, t_IO_final, t_all_frame, t_end_loop, t_init

#-----------------------------------------
PSF = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'files/adk4AKE.psf')))
filenames = os.listdir(os.path.join(os.getcwd(),'files'))
longXTC = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'newtraj.xtc')))
HDF5 = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'newtraj.h5')))
HDF5_1 = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'files/newtraj{}.h5'.format(j))))

if rank == 0:
   copyfile(HDF5, HDF5_1)
   u = mda.Universe(PSF, longXTC)
   print(len(u.trajectory))

MPI.COMM_WORLD.Barrier()

start1 = time.time()
#----------------------------------------------------------------------
u = mda.Universe(PSF, longXTC)
mobile = u.select_atoms("(resid 1:29 or resid 60:121 or resid 160:214) and name CA")
index = mobile.indices
ref0 = mobile
xref0 = ref0.positions-ref0.center_of_mass()
bsize = int(np.ceil(mobile.universe.trajectory.n_frames / float(size)))

f = h5py.File(HDF5, 'r', driver='mpio', comm=MPI.COMM_WORLD)
f.atomic = False
dset = f.require_dataset('pos',(2512208,len(mobile.atoms)*3),dtype='f')

# Create each segment for each process
#for j in range(1,6): # changing files (5 files per block size)
start2 = time.time()

frames_seg = np.zeros([size,2], dtype=int)
for iblock in range(size):
    frames_seg[iblock, :] = iblock*bsize, (iblock+1)*bsize

d = dict([key, frames_seg[key]] for key in range(size))

start, stop = d[rank][0], d[rank][1]
if stop>2512200: stop=2512200

# Block-RMSD in Parallel
start3 = time.time()
out = block_rmsd(dset, xref0, start=start, stop=stop, step=1)

start4 = time.time()
f.close()

start5 = time.time()
if rank == 0:
   data1 = np.zeros([size*bsize,2], dtype=float)
else:
   data1 = None

comm.Gather(out[0], data1, root=0)

start6 = time.time()
if rank == 0:
   data = np.zeros([size,5], dtype=float)
else:
   data = None

comm.Gather(np.array(out[1:], dtype=float), data, root=0)

start7 = time.time()

if rank == 0 and int(j) == 1:
   res = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'RMSD_{}.csv'.format(size))))
   np.save(res, data1)

#print('Cost Calculation')
init_time = start2-start1
comm_time1 = start3-start2
comm_time2 = start6-start5
comm_time3 = start7-start6
comp_time = start4-start3
close_time = start5-start4
tot_time = comp_time+comm_time2

tot_time = comm.gather(tot_time, root=0)
init_time = comm.gather(init_time, root=0)
comm_time1 = comm.gather(comm_time1, root=0)
comm_time2 = comm.gather(comm_time2, root=0)
comm_time3 = comm.gather(comm_time3, root=0)
comp_time = comm.gather(comp_time, root=0)
close_time = comm.gather(close_time, root=0)

if rank == 0:
   tot_time = tot_time
   init_time = init_time
   comm_time1 = comm_time1
   comm_time2 = comm_time2
   comm_time3 = comm_time3
   comp_time = comp_time
   close_time = close_time
else:
   tot_time = None
   init_time = None
   comm_time1 = None
   comm_time2 = None
   comm_time3 = None
   comp_time = None
   close_time = None

# Storing the data
if rank == 0:
   os.remove('files/newtraj{}.h5'.format(j))
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
        file2.write("{} {} {}\n".format(size, j, close_time))

MPI.Finalize
