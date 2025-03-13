import os
import math
import numpy as np
from func import *
import argparse
from matplotlib import pyplot as plt

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',      default=None, type=str)
parser.add_argument('--snapnum',  default=None, type=int)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap splashback features per accretion rate cuts <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Load halos data
# -----------------------------------------------------------------------------------------
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fname = os.listdir(halos_dir)[0]
print(halos_fname)
data = np.load(os.path.join(halos_dir, halos_fname), allow_pickle=True).item()

accretion_rate = data['accretion_rate'] # [dimless]
del data

# Plot the histogram of accretion rates
# -----------------------------------------------------------------------------------------
plt.figure(figsize=(8, 6))
# Calculate the number of halos in each bins
bin_min, bin_max, bin_width = 0, 6, 1
bin_edges = np.arange(bin_min, bin_max + bin_width, bin_width)
bin_centers = bin_edges[:-1] + 0.5 * bin_width
# Plot the histogram
plt.hist(accretion_rate, bins=100, alpha=0.7, color='blue')
plt.xlim(0, 5)
if 'MTNG' in args.sim:
    plot_name = 'MTNG'
elif 'TNG300' in args.sim:
    plot_name = 'TNG300'
plt.savefig(f'result/{plot_name}.png')

# Calculate the number of halos in each bins
hist, _ = np.histogram(accretion_rate, bins=bin_edges)
print(hist)