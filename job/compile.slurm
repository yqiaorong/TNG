#!/bin/bash
#SBATCH --partition=sched_mit_mki_r8  # [sched_mit_mki_r8 / sched_mit_mvogelsb]
#SBATCH --cpus-per-task=4
#SBATCH -o ./output/compile-%j.out
#SBATCH -e ./output/compile-%j.err
#SBATCH --open-mode=append 
#SBATCH --mail-type=BEGIN,END                
#SBATCH --mail-user=s_qyu@mit.edu  

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# Running in directory: "$(pwd)

python3.11 code/DMhalo_MTNG/2_compile_chunk_profile2.py

echo "# Ending at: "$(date)