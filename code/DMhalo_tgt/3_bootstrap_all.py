import os

sim = 'DM-Arepo/MTNG-L500-4320-A/'

snaps = [129, 151, 179, 214, 237, 264]
bin_ends = [4, 4.5, 4.5, 5, 5, 5.5]

for i in range(len(snaps)):

    os.system(f'python3.11 code/DMhalo_tgt/bootstrap_tgt.py --sim {sim} --Nboots 128 '+
              f'--bin_start 3.5 --bin_end {bin_ends[i]} --snapnum {snaps[i]}')