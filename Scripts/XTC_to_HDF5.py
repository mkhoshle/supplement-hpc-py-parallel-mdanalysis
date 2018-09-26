#!/usr/bin/env python
from __future__ import print_function, division
import sys
import numpy as np
import MDAnalysis as mda
import time
from shutil import copyfile
import glob, os
from MDAnalysis import Writer
import mpi4py
from mpi4py import MPI
import h5py 

MPI.Init

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

topdir = os.curdir
PSF = os.path.join(topdir, "files", "adk4AKE.psf")
XTC = os.path.join(topdir, "files", "newtraj.xtc")
HDF5 = os.path.join(topdir, "files",'newtraj.h5')

f = h5py.File(HDF5, 'w', driver='mpio', comm=comm) 
f.atomic = False

u = mda.Universe(PSF, XTC)
mobile = u.select_atoms("(resid 1:29 or resid 60:121 or resid 160:214) and name CA")
topology, trajectory = mobile.universe.filename, mobile.universe.trajectory.filename
index = mobile.indices
bsize = int(np.ceil(u.trajectory.n_frames/float(size)))
frames_seg = np.zeros([size,2], dtype=int)
dfset = f.create_dataset('pos',(bsize*size,len(mobile.atoms)*3),dtype='f')

def get_traj_h5(dfset, index, topology, trajectory, start=None, stop=None, step=None):
    clone = mda.Universe(topology, trajectory)
    g = clone.atoms[index]

    print("block_rmsd", rank, start, stop, step)
    
    frame_i = start
    for iframe, ts in enumerate(clone.trajectory[start:stop:step]):
        dfset[frame_i,:] = np.reshape(g.positions,len(g.atoms)*3)
        frame_i = frame_i+1

    print('rank {} finished working'.format(rank))

    return None

for iblock in range(size):
    frames_seg[iblock, :] = iblock*bsize, (iblock+1)*bsize

d = dict([key, frames_seg[key]] for key in range(size))

start, stop = d[rank][0], d[rank][1]
step = 1
if stop>len(u.trajectory): stop=len(u.trajectory)

out = get_traj_h5(dfset, index, topology, trajectory, start=start, stop=stop, step=step)

comm.Barrier()
f.close()
 
MPI.Finalize
