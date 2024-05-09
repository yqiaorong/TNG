import os
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,type=int)
parser.add_argument('--res',default=1250,type=int)
parser.add_argument('--snap', default=99, type=int)
args = parser.parse_args()

bins = [3.5, 4, 4.5, 5]

# Single stack
for i in range(len(bins)-1):
    os.system(f'python3 code/DMhalo/stacked_density_profiles.py --boxsize {args.boxsize} --res {args.res} --snap {args.snap} --bin_start {bins[i]} --bin_end {bins[i+1]}')

# All stacks together
os.system(f'python3 code/DMhalo/stacked_density_profiles.py --boxsize {args.boxsize} --res {args.res} --snap {args.snap} --bin_start {bins[0]} --bin_end {bins[-1]}')