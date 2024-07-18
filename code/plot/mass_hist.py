import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,type=int)
parser.add_argument('--res',default=1250,type=int)
parser.add_argument('--snapnum',default=99,type=int)
parser.add_argument('--bin_start',default=0,type=float) # 10^{10+x} MSun/h
parser.add_argument('--bin_end',default=6,type=float) # 10^{10+x} MSun/h
args = parser.parse_args()

print('')
print(f'>>> DM halo mass histogram <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Specify the snapshot
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(args.boxsize,args.res)
snapnum = args.snapnum

# Load Halos from groupcat
Group_M_Mean200 = il.groupcat.loadHalos(basePath, snapnum, fields='Group_M_Mean200')


# load hubble param and scale factor
with h5py.File(il.snapshot.snapPath(basePath, snapnum), 'r') as f:
    header = dict(f['Header'].attrs.items())
    scale_factor = header['Time']
    z = np.round( 1 / scale_factor - 1, 3)
    
        

# Save directory
save_dir = f'result/DM halos mass histogram/sim_{args.boxsize}_{args.res}'
if os.path.isdir(save_dir) == False:
    os.makedirs(save_dir)
    
# Histogram of halos mass
plt.figure()
bin_start, bin_end = args.bin_start, args.bin_end
num_bin = int((bin_end-bin_start)*10)
hist_values = plt.hist(Group_M_Mean200, bins=np.logspace(bin_start, bin_end, num_bin))[0]
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