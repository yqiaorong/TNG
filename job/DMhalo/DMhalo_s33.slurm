#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH -t 4320            # Runtime in minutes
#SBATCH -p serial_requeue  
#SBATCH --mem-per-cpu=4000
#SBATCH --open-mode=append # Append when writing files
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/DMhalo/TNG300DM_snap33-%A.out
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/DMhalo/TNG300DM_snap33-%A.err
#SBATCH --mail-type=BEGIN,END    
#SBATCH --mail-user=s_qyu@mit.edu

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# Running in directory: "$(pwd)
    
echo "# Start python processing at: "$(date)
mpirun -np $SLURM_NTASKS python3 code/DMhalo_para/1_subset.py --snapnum 33 --DM _DM --bin_start 1 --bin_end 4.5
echo "# End python processing at: "$(date)

echo "# Ending at: "$(date)