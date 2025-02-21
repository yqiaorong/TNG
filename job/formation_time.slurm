#!/bin/bash
#SBATCH --ntasks-per-node=11
#SBATCH -t 4320            
#SBATCH -p serial_requeue  
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-10
#SBATCH --open-mode=append 
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/formation_time/TNG300-%A-%a.out
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/formation_time/TNG300-%A-%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=s_qyu@mit.edu

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# Running in directory: "$(pwd)

snaplist=( $(cat job/boots_TNG300/snap_list.txt) )
snap=${snaplist[$SLURM_ARRAY_TASK_ID]}

echo "# Array ID: " $SLURM_ARRAY_TASK_ID
    
echo "# Start python processing $SLURM_ARRAY_TASK_ID at: "$(date)
python3 /n/home01/sqyu/AstroLab/TNG/code/DMhalo_TNG300/5_add_formation_time.py --snapnum $snap --sim TNG300/sim_205_1250_DM/
echo "# End python processing $SLURM_ARRAY_TASK_ID at: "$(date)

echo "# Ending at: "$(date)
