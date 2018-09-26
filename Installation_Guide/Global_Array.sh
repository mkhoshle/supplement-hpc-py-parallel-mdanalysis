# ---------Notes---------------------------
# There are many configurations of the GA build.
# By default, it will use “MPI-TS” (--with-mpi-ts) during configure.
# This is the most compatible build for GA, but not the fastest.
# The fastest is MPI-PR (--with-mpi-pr).
# There is also the “MPI-MT” (--with-mpi-mt) and “MPI-PT” (--with-mpi-pt) configure options.  
# Those perform better than TS, in general, but still lag behind PR.  
# In an infiniband cluster, --with-openib can be used during configure.

## REQUIREMENTS: Cython, Numpy, Openmpi, MPI4py, Gfortran ##
###### COMET #######
[1] Set up the MPI environment:

conda create -n ga4py python=2.7

[2] Steps for the conda based install: Scipy Numpy

conda install scipy numpy ipython

[3] Steps for installing Openmpi:
tar xf openmpi-1.10.7.tar.gz
cd openmpi-1.10.7
mkdir /home/mkhoshle/openmpi-1.10.7-install
./configure --prefix=/home/mkhoshle/openmpi-1.10.7-install --with-devel-headers --enable-binaries
make -j 4
make install

export LD_LIBRARY_PATH="/home/mkhoshle/openmpi-1.10.7-install/lib:$LD_LIBRARY_PATH"
export PATH="/home/mkhoshle/openmpi-1.10.7-install/bin:$PATH"

[4] Install MPI4py:

pip install mpi4py

[5] Steps for installing Global Array:

tar xf ga-5.6.1.tar.gz
cd ga-5.6.1
mkdir /home/mkhoshle/ga-5.6.1-install
./configure --prefix=/home/mkhoshle/ga-5.6.1-install --with-openib --disable-static --enable-shared --enable-i8
make -j 4
make install

export LD_LIBRARY_PATH="/home/mkhoshle/ga-5.6.1-install/lib:$LD_LIBRARY_PATH"
export PATH="/home/mkhoshle/ga-5.6.1-install/bin:$PATH"

[6] Install Cython:

pip install cython

[7] Install Ga4py
git clone https://github.com/GlobalArrays/ga4py
cd ga4py
python setup.py build 
python setup.py install --prefix=/home/mkhoshle/miniconda2/envs/ga4py

[8] Set the path for precompiled gfortran

tar xf gcc-4.9.4.tar.xz
export LD_LIBRARY_PATH="/home/mkhoshle/gcc-4.9.4/lib64:$LD_LIBRARY_PATH"
export PATH="/home/mkhoshle/gcc-4.9.4/bin:$PATH"

export PYTHONPATH="/home/mkhoshle/miniconda2/envs/ga4py/lib/python2.7/site-packages"
