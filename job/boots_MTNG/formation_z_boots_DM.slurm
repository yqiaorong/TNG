#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH -t 4320            
#SBATCH -p serial_requeue  
#SBATCH --mem-per-cpu=4000
#SBATCH --open-mode=append 
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/boots_formation_z/MTNG-DM-%A-%a.out
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/boots_formation_z/MTNG-DM-%A-%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=s_qyu@mit.edu

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# Running in directory: "$(pwd)
    
python3 /n/home01/sqyu/AstroLab/TNG/code/bootstrap/bootstrap_z.py --snapnum 264 --sim MTNG/DM-Arepo/MTNG-L500-4320-A/

echo "# Ending at: "$(date)
