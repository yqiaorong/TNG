#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2  
#SBATCH -t 4320              # Runtime in minutes
#SBATCH -p serial_requeue   
#SBATCH --mem-per-cpu=4000
#SBATCH --open-mode=append   # Append when writing files
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/tree_snap_mass_idx_%j.out 
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/tree_snap_mass_idx_%j.err 
#SBATCH --mail-type=BEGIN,END    
#SBATCH --mail-user=s_qyu@mit.edu

python3 code/tree/get_halo_idx_table_at_snap_X.py