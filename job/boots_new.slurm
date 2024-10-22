#!/bin/bash
#SBATCH --ntasks-per-node=11
#SBATCH -t 4320            
#SBATCH -p serial_requeue  
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-10
#SBATCH --open-mode=append 
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/boots/TNG300-Hydro-%A-%a.out
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/boots/TNG300-Hydro-%A-%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=s_qyu@mit.edu

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# Running in directory: "$(pwd)

snaplist=( $(cat job/snap_list.txt) )
snap=${snaplist[$SLURM_ARRAY_TASK_ID]}

bin_end_list=( $(cat job/bin_ends.txt) )
bin_end=${bin_end_list[$SLURM_ARRAY_TASK_ID]}

bin_start_list=( $(cat job/bin_starts.txt) )
bin_start=${bin_start_list[$SLURM_ARRAY_TASK_ID]}

echo "# Array ID: " $SLURM_ARRAY_TASK_ID
    
echo "# Start python processing $SLURM_ARRAY_TASK_ID at: "$(date)
python3 code/DMhalo_tgt/bootstrap_tgt_phys.py --snapnum $snap --bin_start $bin_start --bin_end $bin_end --sim TNG300/sim_205_1250_Hydro/
echo "# End python processing $SLURM_ARRAY_TASK_ID at: "$(date)

echo "# Ending at: "$(date)