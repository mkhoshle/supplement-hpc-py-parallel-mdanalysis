# Supplementary Material for _Parallel Performance of Molecular Dynamics Trajectory Analysis_

[![DOI](https://zenodo.org/badge/146852628.svg)](https://zenodo.org/badge/latestdoi/146852628)

Supplementary material (scripts and documentation) for the paper _Parallel Performance of Molecular Dynamics Trajectory Analysis_, available as preprint [arXiv:1907.00097](https://arxiv.org/abs/1907.00097).

For paper source files and data for figures see repository https://github.com/hpcanalytics/paper-hpc-py-parallel-mdanalysis. Each release of this repository is archived in the Zenodo repository under DOI [10.5281/zenodo.3351616](https://doi.org/10.5281/zenodo.3351616).

# Contents and License

- benchmark scripts
- notes on setting up the computational environment on the [XSEDE](https://www.xsede.org/) resources discussed in the paper
- raw data, analysis and plotting code for the figures in the paper

All content is made available under the [GNU General Public
License v3](https://www.gnu.org/licenses/gpl.html) (code) or the
[CC-BY-SA
4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode) (text
and documentation).

# Requirements
* Python 2.7 or higher or Python 3.1 or higher 
* GCC 4.9.4  
* OpenMPI 1.10.7
* MPI4py 3.0.0
* PHDF5 1.10.1
* H5py 2.7.1
* MDAnalysis 0.17.0 or higher

# Some Comments on Installation

Different packages and libraries are used in the present study in order to achieve the desired performance.
From the libraries used in the present study (given in Table 2 of the manuscript), _MPI4py_, _H5py_, _MDAnalysis_ were the easiest to build. 
Users are able to get these packages installed through the _Conda_ package manager [Conda](https://conda.io/docs/).

_OpenMPI_, and _GCC_, can be easily configured and installed.
The configure script supports many different command line options, but the support for these libraries is very strong, they are widely used and the excellent documentation and discussion mailing lists provide users with great resources for consult, troubleshooting and tracking issues.

 _PHDF5_ is harder to install especially because the libraries are built from source and variable configure options are required to support specific network interconnects and back-end run-time environments.

We performed our benchmark on several HPC resources and detailed instruction for installing all these packages on different resources are provided.

# Getting things running
To have things running efficiently you may need a lot of additional setup depending on the resource. Take a look at  `Installation_Guide/h5py-phdf5.sh` for examples which explains how we have setup the computational evironment in our study. 

We have provided all the steps required for building the packages and running sample test cases on all the XSEDE resources used in our study. For the cases which the instructions are given on one resource only, the steps have been the same on other resources as well.

# Scripts
All the scripts for running `RMSD` benchamarks can be found at `Scripts` directory.


