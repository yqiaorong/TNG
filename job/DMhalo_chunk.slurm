#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1                   # adjust for hardware (roughly #cores/#gpus)
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=sched_mit_mki_r8
#SBATCH --nodelist=node2001                 # Specify the node
#SBATCH -o ./output/DMhalo/halo-DM214_chunk-%A.out
#SBATCH -e ./output/DMhalo/halo-DM214_chunk-%A.err
#SBATCH --mail-type=BEGIN,END                
#SBATCH --mail-user=s_qyu@mit.edu  

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# GPU spec: "$(nvidia-smi --query-gpu=gpu_name,gpu_bus_id,gpu_uuid --format=csv,noheader)
echo "# Running in directory: "$(pwd)

echo "# Start python processing at: "$(date)
python3.11 code/DMhalo_MTNG/1_serial.py --sim DM-Arepo/MTNG-L500-4320-A --snapnum 214 --bin_start 3 --start_chunk_idx 159 --start_chunk_idx 0
echo "# End python processing at: "$(date)

echo "# Ending at: "$(date)