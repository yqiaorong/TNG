#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2  
#SBATCH -t 4320               # Runtime in minutes
#SBATCH -p serial_requeue 
#SBATCH --mem-per-cpu=4000
#SBATCH --open-mode=append    # Append when writing files
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/compile_single_%j.out 
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/compile_single_%j.err 
#SBATCH --mail-type=BEGIN,END    
#SBATCH --mail-user=s_qyu@mit.edu

echo "# Start python processing at: "$(date)
python3 code/DMhalo/compile_chunk_profile.py --snapnum 99 --DM _DM --bin_start 1 --bin_end 4.5
echo "# End python processing at: "$(date)