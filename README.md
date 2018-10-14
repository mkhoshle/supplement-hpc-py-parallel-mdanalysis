# Supplementary Material

Supplementary material (scripts and documentation) for the paper _Guidelines on Obtaining Ideal Performance for Parallel Analysis of Molecular Dynamics Trajectories with Python on HPC Resources_ (see repository https://github.com/hpcanalytics/paper-hpc-py-parallel-mdanalysis).

- benchmark scripts
- notes on setting up the computational environment on the [XSEDE](https://www.xsede.org/) resources discussed in the paper

All content is made available under the [GNU General Public
License v3](https://www.gnu.org/licenses/gpl.html) (code) or the
[CC-BY-SA
4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode) (text
and documentation).

# Requirements
* Python 2.7 or higher or Python 3.1 or higher 
* GCC 4.9.4  
* OpenMPI 1.10.7
* Global Array 5.6.1
* MPI4py 3.0.0
* GA4py 1.0
* PHDF5 1.10.1
* H5py 2.7.1
* MDAnalysis 0.17.0 or higher

# Some Comments on Installation

Different packages and libraries are used in the present study in order to achieve the desired performance.
From the libraries used in the present study (given in Table 2 of the manuscript), _MPI4py_, _H5py_, _MDAnalysis_ were the easiest to build. 
Users are able to get these packages installed through the _Conda_ package manager [Conda](https://conda.io/docs/).

_OpenMPI_, and _GCC_, can be easily configured and installed.
The configure script supports many different command line options, but the support for these libraries is very strong, they are widely used and the excellent documentation and discussion mailing lists provide users with great resources for consult, troubleshooting and tracking issues.

_Global Array_, and _PHDF5_ are harder to install especially because the libraries are built from source and variable configure options are required to support specific network interconnects and back-end run-time environments.
_GA4py_ has only one release, and does not provide users with extensive documentation.

We performed our benchmark on several HPC resources and detailed instruction for installing all these packages on different resources are provided.

# Getting things running
To have things running efficiently you may need a lot of additional setup depending on the resource. Take a look at `Installation_Guide/Global_Array.sh` and `Installation_Guide/h5py-phdf5.sh` for examples which explains how we have setup the computational evironment in our study. 

We have provided all the steps required for building the packages and running sample test cases on all the XSEDE resources used in our study. For the cases which the instructions are given on one resource only, the steps have been the same on other resources as well.

# Scripts
All the scripts for running `RMSD` benchamarks can be found at `Scripts` directory.


