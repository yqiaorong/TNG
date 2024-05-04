import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--snap',default=99,type=int)
args = parser.parse_args()

print('')
print(f'>>> DM halo mass histogram <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Specify the snapshot
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L205n1250TNG/output'
snapNum = args.snap

# Load Halos from groupcat
Group_M_Mean200 = il.groupcat.loadHalos(basePath, snapNum, fields='Group_M_Mean200')

# Save directory
save_dir = 'result/DM halos mass histogram'
if os.path.isdir(save_dir) == False:
    os.makedirs(save_dir)
    
# Histogram of halos mass
plt.figure()
hist_values, _, _ = plt.hist(Group_M_Mean200, bins=np.logspace(0.01, 6, 50))
plt.xlabel('Mass [10^10 MSun / h]')
plt.ylabel('Frequency')
plt.yscale('log')
plt.xscale('log')
plt.title(f"DM halos' mass histogram at snap {snapNum}")
plt.savefig(f'result/DM halos mass histogram/snap_{snapNum}')

# Print numbers in each bin
sums = []
for i in range(0, len(hist_values), 5):
    sum_of_five = sum(hist_values[i:i+5])
    sums.append(sum_of_five)
print(sums)