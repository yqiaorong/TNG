#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH -t 4320          
#SBATCH -p serial_requeue
#SBATCH --mem-per-cpu=4000
#SBATCH --open-mode=append 
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/bootstrap_all_%j.out 
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/bootstrap_all_%j.err 
#SBATCH --mail-type=BEGIN,END                
#SBATCH --mail-user=s_qyu@mit.edu  

mpirun -np $SLURM_NTASKS python3 code/DMhalo_para/3_bootstrap_all.py 