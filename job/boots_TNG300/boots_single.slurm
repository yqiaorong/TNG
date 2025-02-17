#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH -t 4320            
#SBATCH -p serial_requeue  
#SBATCH --mem-per-cpu=4000
#SBATCH --open-mode=append 
#SBATCH -o /n/home01/sqyu/AstroLab/TNG/output/boots_mass/TNG300-DM-%A-%a.out
#SBATCH -e /n/home01/sqyu/AstroLab/TNG/output/boots_mass/TNG300-DM-%A-%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=s_qyu@mit.edu

echo "# Starting at: "$(date)
echo "# Running on machine: "$(hostname)
echo "# Running in directory: "$(pwd)

python3 /n/home01/sqyu/AstroLab/TNG/code/bootstrap/bootstrap_mass.py --snapnum 25 --sim TNG300/sim_205_1250_DM/

echo "# Ending at: "$(date)
