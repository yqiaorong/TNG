"""This scripts plot the distribution of bootstrapped splashback features. 
The target is to verify if they follow normal distribution."""

import os
import numpy as np
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')

Nboots = 128
mass_cut_idx = 0 # 10^11 Msun/h
feat_idx = 2 # width dimensionless
root_dir = 'result/bootstrap_stats/'

save_dir = f'result/hist/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
# Load data directory
TNG300_dir = f'{root_dir}/sim_205_1250/Nboots_{Nboots}/'
TNG300_list = os.listdir(TNG300_dir)

TNG300_cmap = plt.get_cmap('winter', len(TNG300_list))
TNG300_snaps = [99, 78, 67, 50, 40, 33, 21, 17, 13, 8]
TNG300_mass_cuts = [10**11, 10**11.5, 10**12, 10**12.5, 10**13, 10**13.5]

# Plot 1
fig, axs = plt.subplots(1, len(TNG300_snaps), figsize=(10, 2), dpi=500)
for isnap, snap in enumerate(TNG300_snaps):
    
    data = np.load(TNG300_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data = data['full_results'][mass_cut_idx, feat_idx, :]
    axs[isnap].hist(data)
    
plt.savefig(f'{save_dir}/distri_mass{mass_cut_idx}_feat{feat_idx}')
plt.close()

# Plot 2
fig, axs = plt.subplots(1, 1, dpi=500)
for isnap, snap in enumerate(TNG300_snaps):
    
    data = np.load(TNG300_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data = data['full_results'][mass_cut_idx, feat_idx, :]
    axs.hist(data)
    
plt.savefig(f'{save_dir}/distri_mass{mass_cut_idx}_feat{feat_idx}_full')
plt.close()

# Plot 3
for isnap, snap in enumerate(TNG300_snaps):
    
    fig, axs = plt.subplots(1, 1, dpi=500)
    
    data = np.load(TNG300_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data = data['full_results'][mass_cut_idx, feat_idx, :]
    axs.hist(data)
    
    plt.savefig(f'{save_dir}/distri_snap{snap}_mass{mass_cut_idx}_feat{feat_idx}')
    plt.close()