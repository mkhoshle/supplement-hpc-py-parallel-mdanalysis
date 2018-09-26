###### Bridge ######
[1] Set up the MPI environment:

conda create -n mpio python=3
source activate mpio
module purge
module list
module load phdf5/1.8.16_intel

[2] Install ipython, mpi4py, mdanalysis:

pip install ipython
pip install mpi4py
pip install --upgarde mdanalysis

[3] Install h5py:

export CC="/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpicc"
export HDF5_MPI="ON"
pip install --no-binary=h5py h5py

###### COMET #######
[1] Set up the MPI environment:

module purge
module load gnutools
module load gnu
export MODULEPATH=/home/mkhoshle/modulefiles:$MODULEPATH
module load openmpi_ib/1.10.7 # Instruction for Openmpi installation are given in Global_Array.sh

[2] Configure:

tar -xf hdf5-1.10.1.tar.gz
cd hdf5-1.10.1/
export CC=mpicc
export FC=mpif90
export CXX=mpicxx
export CFLAGS="-O2 -mavx"
export FCFLAGS="-O2 -mavx"
export CXXFLAGS="-O2 -mavx"
./configure --prefix=/home/mkhoshle/hdf5/v1.10.1 --enable-parallel

[3] Build and install
make
make install

[4] Setup the module files in the home directory.

(a) Create the directory:

/home/mkhoshle/modulefiles/hdf5

(b) Create the module and version files:

1.10.1
.version.1.10.1

#-------------------------------------------------------------------
module show hdf5/1.10.1
-------------------------------------------------------------------
/home/mkhoshle/modulefiles/hdf5/1.10.1:

module-whatis	 hdf5 
module-whatis	 Version: 1.10.1 
module-whatis	 Description: hdf5 
module-whatis	 Compiler: gnu 
module-whatis	 MPI Flavors: openmpi_ib 
setenv		 HDF5HOME /home/mkhoshle/hdf5/v1.10.1 
setenv		 HDF5_DIR /home/mkhoshle/hdf5/v1.10.1 
setenv		 HDF5_ROOT /home/mkhoshle/hdf5/v1.10.1 
setenv		 PHDF5_VERSION 1.10.1 
setenv		 PHDF5_HOME /home/mkhoshle/hdf5/v1.10.1 
setenv		 PHDF5_ROOT /home/mkhoshle/hdf5/v1.10.1 
prepend-path	 PATH /home/mkhoshle/hdf5/v1.10.1/bin 
prepend-path	 LD_LIBRARY_PATH /home/mkhoshle/hdf5/v1.10.1/lib 
prepend-path	 LIBPATH /home/mkhoshle/hdf5/v1.10.1/lib 
#----------------------------------------------------------------------------------------------------

[1] Make sure you do any of your conda installs on a compute node (Otherwise you will get errors due to the thread limit issue on the login node).

 srun --nodes=1 --ntasks-per-node=24 -p debug -t 00:30:00 --pty --wait 0 /bin/bash

[2] Once on the compute node:

module purge
module load gnutools
module load gnu
export MODULEPATH=/home/mkhoshle/modulefiles:$MODULEPATH
module load openmpi_ib/1.10.7
module load hdf5/1.10.1

[3] Steps for the conda based install:

conda create -n mpio python=3
source activate mpio
pip install ipython
pip install mpi4py
pip install --upgarde MDAnalysis
export CC="/home/mkhoshle/openmpi-1.10.7-install/bin/mpicc"
export HDF5_MPI="ON"
pip install --no-binary=h5py h5py

#For Running a sample case:
[1] Request an interactive node:

 srun --nodes=1 --ntasks-per-node=24 -p debug -t 00:30:00 --pty --wait 0 /bin/bash

[2] Set up the environment:

module purge
module load gnutools
module load gnu
export MODULEPATH=/home/mkhoshle/modulefiles:$MODULEPATH
module load openmpi_ib/1.10.7
module load hdf5/1.10.1
source activate mpio

[3] Run the code:

mpirun -x PATH -x LD_LIBRARY_PATH --mca mpi_warn_on_fork 0 -n 4 python test.py

[4] Check the output

h5dump parallel_test.hdf5 #Gives sample output. 


