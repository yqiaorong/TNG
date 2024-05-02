#!/bin/bash

# Job Name and Files
#SBATCH -J DM_halo

#Output and error
#SBATCH -o ./output/DM_halo_%j.out
#SBATCH -e ./output/DM_halo_%j.err
#Initial working directory
#SBATCH -D ./

#Partition & time limit
#SBATCH --partition=sched_mit_mvogelsb
#SBATCH --time=96:00:00

#Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#SBATCH --constraint=centos7
#SBATCH --mem-per-cpu=4000 # 4GB of memory per CPU

#SBATCH --exclusive
#SBATCH --export=ALL

#SBATCH --mail-type=ALL
#SBATCH --mail-user=s_qyu@mit.edu
module purge

module load gcc/9.3.0
module load openmpi/4.0.5
module load fftw/3.3.9
module load hdf5_18/1.8.12
module load gsl/2.5

python3 code/DMhalo/1_all_halos.py