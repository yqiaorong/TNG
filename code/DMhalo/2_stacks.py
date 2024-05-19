import os
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,type=int)
parser.add_argument('--res',default=1250,type=int)
args = parser.parse_args()

snaps = [25, 33, 40, 50, 67, 78, 99]

# ### Bin 3 - 3.5
# # Single stack
# for s in snaps:
#     os.system(f'python3 code/DMhalo/stacked_density_profiles.py'+
#               f' --boxsize {args.boxsize} --res {args.res} --snap {s}'+
#               f' --bin_start 3 --bin_end 3.5'+
#               f' --root_dir DMhalo_density_profiles_old')
# # Profile time evolution
# os.system(f'python3 code/DMhalo/profiles_time_evolution.py'+
#           f' --boxsize {args.boxsize} --res {args.res}'+
#           f'--bin_start 3 --bin_end 3.5'+
#           f' --root_dir DMhalo_density_profiles_old')

### Bin 3.5 - 4
# Single stack
for s in snaps[1:]:
    os.system(f'python3 code/DMhalo/stacked_density_profiles.py'+
              f' --boxsize {args.boxsize} --res {args.res} --snap {s}'+
              f' --bin_start 3.5 --bin_end 4'+
              f' --root_dir DMhalo_density_profiles_old')
# Profile time evolution
os.system(f'python3 code/DMhalo/profiles_time_evolution.py'+
          f' --boxsize {args.boxsize} --res {args.res}'+
          f' --bin_start 3.5 --bin_end 4'+
          f' --root_dir DMhalo_density_profiles_old')