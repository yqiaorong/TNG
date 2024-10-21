#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2  # n core should line up with the number of snaps here
#SBATCH -t 4320 # Runtime in minutes
#SBATCH -p serial_requeue    # Partition to submit to
#SBATCH --mem-per-cpu=4000
#SBATCH --open-mode=append   # Append when writing files
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/tree_mass_table_%j.out 
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/tree_mass_table_%j.err 
#SBATCH --mail-type=BEGIN,END    
#SBATCH --mail-user=s_qyu@mit.edu

python3 code/tree/mass_table.py --DM _DM