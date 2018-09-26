#!/usr/bin/env python
from __future__ import print_function, division
import sys
import numpy as np
import MDAnalysis as mda
import time
import glob, os
from MDAnalysis import Writer
import mpi4py
from mpi4py import MPI

MPI.Init

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

PSF = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'files/adk4AKE.psf')))
XTC = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(),'newtraj.xtc')))

u = mda.Universe(PSF, XTC)

bsize = int(np.ceil(u.trajectory.n_frames/float(size)))
frames_seg = np.zeros([size,2], dtype=int)
for iblock in range(size):
    frames_seg[iblock, :] = iblock*bsize, (iblock+1)*bsize

d = dict([key, frames_seg[key]] for key in range(size))

start, stop = d[rank][0], d[rank][1]
step = 1

XTC_seg = os.path.abspath(os.path.normpath(os.path.join(os.getcwd(), "traj_{}".format(size),"newtraj_{}.xtc".format(rank))))

start1 = time.time()
with mda.Writer(XTC_seg) as out:
    for iframe, ts in enumerate(u.trajectory[start:stop:step]):
        out.write(ts)

start2 = time.time()

print("{} for rank {}".format(start2-start1,rank))

MPI.Finalize
