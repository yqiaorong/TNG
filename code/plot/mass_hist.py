import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',    default=None, type=str)
parser.add_argument('--snapnum',default=264,  type=int) 
# Available snapshots: [264, 237, 214, 179, 151, 129, 094, 080, 069, 051]
parser.add_argument('--bin_start',default=1.0,type=float) # 10^{10+x} MSun/h
parser.add_argument('--bin_end',  default=6.0,type=float) # 10^{10+x} MSun/h
args = parser.parse_args()

print('')
print(">>> DM halo mass histogram <<<")
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Specify the snapshot
data_path = '/virgotng/mpa/MTNG/'
basePath = data_path + args.sim
snapnum = args.snapnum

# Load Halos from groupcat
Group_M_Mean200 = il.groupcat.loadHalos(basePath, snapnum, fields='Group_M_Mean200')


# load hubble param and scale factor
with h5py.File(il.snapshot.snapPath(basePath, snapnum), 'r') as f:
    header = dict(f['Header'].attrs.items())
    scale_factor = header['Time']
    z = np.round( 1 / scale_factor - 1, 3)
    
        

# Save directory
save_dir = f"result/DM_halos_mass_histogram/{args.sim}"
if os.path.isdir(save_dir) == False:
    os.makedirs(save_dir)
    
# Histogram of halos mass
plt.figure()
bin_start, bin_end = args.bin_start, args.bin_end
num_bin = int((bin_end-bin_start)*10)
bins = np.logspace(bin_start, bin_end, num_bin)
hist_values = plt.hist(Group_M_Mean200, bins=bins)[0]
plt.xlabel('Mass [10^10 MSun / h]')
plt.ylabel('Frequency')
plt.yscale('log')
plt.xscale('log')
plt.title(f"DM halos' mass histogram at z = {z}")
plt.savefig(os.path.join(save_dir, f'snap_{snapnum}'))

# Print numbers in each bin
sums = []
for i in range(0, num_bin, 5):
    sum_of_five = sum(hist_values[i:i+5])
    sums.append(sum_of_five)
print(sums)
