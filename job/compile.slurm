#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2  
#SBATCH -t 4320               # Runtime in minutes
#SBATCH -p serial_requeue 
#SBATCH --mem-per-cpu=4000
#SBATCH --open-mode=append    # Append when writing files
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/compile_%j.out 
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/compile_%j.err 

python3 code/DMhalo/2_compile_all.py 