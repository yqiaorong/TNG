import os
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',    default=None, type=str)
parser.add_argument('--bin_start',default=2.0,type=float) # 10^{10+x} MSun/h
parser.add_argument('--bin_end',  default=6.0,type=float) # 10^{10+x} MSun/h
args = parser.parse_args()

snaps = [264, 237, 214, 179, 151, 129, 94, 80, 69, 51]

for snap in snaps:
    os.system(f'python3 code/plot/mass_hist.py --sim {args.sim} --snapnum {snap} --bin_start {args.bin_start} --bin_end {args.bin_end}')
