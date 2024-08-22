#!/bin/bash

# Slurm job options
#SBATCH --account=d114-KFUPM  
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --job-name=Output
#SBATCH --time=08:00:00
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

# This module needs to be loaded in ALL scripts as the first module
module load epcc-job-env

# Now load the "lammps" package
export PATH="$PATH:/work/d114/d114/yichongc/mylammps/src/" 

# Set the number of threads to 1
export OMP_NUM_THREADS=1

# Run the parallel job (input file should be in cwd)
srun -n 512 lmp_archer2 < in.equil
#    Using 512 MPI processes 
#    128 MPI processes per node
