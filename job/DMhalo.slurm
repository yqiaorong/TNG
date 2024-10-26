#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10                  # adjust for hardware (roughly #cores/#gpus)
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-15                        # Adjust number of total job tasks you want to run at once
#SBATCH --partition=sched_mit_mki_r8
#SBATCH --nodelist=node2001                 # Specify the node
#SBATCH -o ./output/DMhalo/halo-DM264-%A-%a.out
#SBATCH -e ./output/DMhalo/halo-DM264-%A-%a.err
#SBATCH --mail-type=BEGIN,END                
#SBATCH --mail-user=s_qyu@mit.edu  

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# GPU spec: "$(nvidia-smi --query-gpu=gpu_name,gpu_bus_id,gpu_uuid --format=csv,noheader)
echo "# Running in directory: "$(pwd)

chunk_list=( $(cat job/DMhalo_chunk.txt) )
chunk_idx=${chunk_list[$SLURM_ARRAY_TASK_ID]}

echo "# Array ID: " $SLURM_ARRAY_TASK_ID

for (( i=0; i<40; i+=4 )); do
    core_start_idx=$((chunk_idx + i))
    
    echo "# Start python processing $core_start_idx at: "$(date)
    python3.11 code/DMhalo/1_serial.py --sim DM-Arepo/MTNG-L500-4320-A --snapnum 264 --bin_start 2 --bin_end 3 --start_chunk_idx $core_start_idx &
    echo "# End python processing $core_start_idx at: "$(date)
done

wait

echo "# Ending at: "$(date)