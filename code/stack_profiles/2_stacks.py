import os
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,type=int)
parser.add_argument('--res',default=1250,type=int)
args = parser.parse_args()

snaps = [17, 21, 25, 33, 40, 50, 67, 78, 99]
bins = [2]

for i in range(len(bins)):
    # Single stack
    for s in snaps:
        os.system(f'python3 code/DMhalo/stacked_density_profiles.py'+
                f' --boxsize {args.boxsize} --res {args.res} --snap {s}'+
                f' --bin_start {bins[i]} --bin_end {bins[i]+0.5}'+
                f' --root_dir DMhalo_density_profiles_old')
    # Profile time evolution
    os.system(f'python3 code/DMhalo/profiles_time_evolution.py'+
            f' --boxsize {args.boxsize} --res {args.res}'+
            f' --bin_start {bins[i]} --bin_end {bins[i]+0.5}'+
            f' --root_dir DMhalo_density_profiles_old')
   
# Sort data
os.system(f'python3 code/DMhalo/Rsp_data_sort.py --root_dir DMhalo_density_profiles_old')