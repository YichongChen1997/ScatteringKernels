# ScatteringKernels
Source code for the project of scattering kernels

There are two types surfaces have been considered, namely plantium and kerogen surfaces.


1. initialisation: contain the C++ code for simulation initialistion. You can specify the channel size, roughness and thickness of the surface. It will generate the data.dat file, which is the input for the equilibration run (in.equil).

2. in.equil : equilibration run, which will generate the dumpe file and restart file.

3. in.meas  : measurement or production run. It will generate the final dump file after the analysis.

4. src      : contain all the post-processing codes for the transport and scattering analysis. 

4. Makefile : compile the C++ codes in a quicker way.

5. equil.slurm/meas.slurm : bash script for submitting jobs in ARCHER2 (HPC).
