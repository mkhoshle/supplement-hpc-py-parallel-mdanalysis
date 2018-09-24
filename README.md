# Supplementary Material

Supplementary material (scripts and documentation) for the paper _Guidelines on Obtaining Ideal Performance for Parallel Analysis of Molecular Dynamics Trajectories with Python on HPC Resources_ (see repository https://github.com/hpcanalytics/paper-hpc-py-parallel-mdanalysis).

- benchmark scripts
- notes on setting up the computational environment on the [XSEDE](https://www.xsede.org/) resources discussed in the paper

All content is made available under the [GNU General Public
License v3](https://www.gnu.org/licenses/gpl.html) (code) or the
[CC-BY-SA
4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode) (text
and documentation).


# Installation

Different packages and libraries are used in the present study in order to achieve the desired performance.
From the libraries given in Table , _MPI4py_, _H5py_, _MDAnalysis_ were the easiest to build. 
Users are able to get these packages installed through the _Conda_ package manager [Conda](https://conda.io/docs/).
_OpenMPI_, and _GCC_, can be easily configured and installed.
The configure script supports many different command line options, but the support for these libraries is very strong, they are widely used and the excellent documentation and discussion mailing lists provide users with great resources for consult, troubleshooting and tracking issues.

_Global Array_, and _PHDF5_ are harder to install especially because the libraries are built from source and variable configure options are required to support specific network interconnects and back-end run-time environments.
_GA4py_ has only one release, and does not provide users with extensive documentation.

We performed our benchmark on several HPC resources and therefore, we had to install all the related packages and tools on all resources.
However, there are always differences in the resources because their set up and architectures differ from each other. 
For example, on SuperMIC although tool installation was done in the same way as Comet and also passed initial testing, the execution did not distribute the processes to all nodes. 
This was due to the fact that our custom OpenMPI installation did not correctly parse the node list offered by SuperMIC to our job. 
Thus, we had to manually pass the node lists to MPIRUN. 
In addition, we found that the loaded modules, along with library path changes, did not propagate to all nodes from our OpenMPI installation. 
OpenMPI's execution engine could not access the correct libraries and was not able to launch the processes correctly. 
Reinstalling OpenMPI with enabling the flag to use the Open Run-Time Environment (ORTE) by default and including the OpenMPI installation to BASHRC allowed for correct execution.
