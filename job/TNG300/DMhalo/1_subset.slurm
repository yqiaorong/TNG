#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH -t 4320         
#SBATCH -p serial_requeue  
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-10
#SBATCH --open-mode=append # Append when writing files
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/DMhalo/TNG300-%A-%a.out
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/DMhalo/TNG300-%A-%a.err
#SBATCH --mail-type=BEGIN,END    
#SBATCH --mail-user=s_qyu@mit.edu

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# Running in directory: "$(pwd)

snaplist=( $(cat job/boots_TNG300/snap_list.txt) )
snap=${snaplist[$SLURM_ARRAY_TASK_ID]}

bin_end_list=( $(cat job/bin_ends.txt) )
bin_end=${bin_end_list[$SLURM_ARRAY_TASK_ID]}

echo "# Array ID: " $SLURM_ARRAY_TASK_ID
    
echo "# Start python processing $SLURM_ARRAY_TASK_ID at: "$(date)
mpirun -np $SLURM_NTASKS python3 code/DMhalo_TNG300/1_subset.py --snapnum $snap --DM DM --bin_start 1 --bin_end $bin_end 
echo "# End python processing $SLURM_ARRAY_TASK_ID at: "$(date)

echo "# Ending at: "$(date)