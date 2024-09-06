import pandas as pd
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import h5py
import illustris_python as il
import numpy as np
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM', default='', type=str)
args = parser.parse_args()



DM = args.DM
boxsize, res = 205, 1250

data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)

load_dir = f'result/DMhalo_mass_table/sim_{boxsize}_{res}{DM}/'

mass_df = pd.read_csv(load_dir+'mass_table.csv', index_col=0)
accretion_df = pd.read_csv(load_dir+'accretion_table.csv', index_col=0)

snap_list = mass_df.index

# Load redshift
redshifts = []
for snap in snap_list:
    with h5py.File(il.snapshot.snapPath(basePath, int(snap[5:])), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        redshifts.append(z)

for isnap in range(len(snap_list)):
    if isnap != 0:
        rate = accretion_df.iloc[isnap, :]
        rate = rate.replace([np.inf, -np.inf], np.nan)
        rate = rate.dropna().to_numpy()
        median = np.median(rate)

        plt.figure()
        hist = plt.hist(rate, label=f'z = {np.round(redshifts[isnap], 2)}')[0]
        plt.plot([median, median], [0, max(hist)], linestyle='dashed', color='red', 
                  label=f'median = {np.round(median, 3)}')
        plt.legend(loc='best')
        plt.savefig(load_dir+f'snap_{snap_list[isnap][5:]}')

# Plot per mass stacks
mass_cuts = np.arange(1, 4.5, 0.5)
num_cuts = int((4.5-1)/0.5)

for isnap in range(len(snap_list)):
    
    if isnap != 0:
        fig, axes = plt.subplots(1, num_cuts, figsize=(15, 5))
        
        mass = mass_df.iloc[isnap, :]
        rate = accretion_df.iloc[isnap, :]
        
        # Drop inf and nan
        rate = rate.replace([np.inf, -np.inf], np.nan)
        valid_mask = ~np.isnan(rate)
        rate, mass = rate[valid_mask], mass[valid_mask]

        for icut, cut in enumerate(mass_cuts[:-1]):
            # Get the mass bin
            mass_mask = (mass >= 10**mass_cuts[icut]) & (mass < 10**mass_cuts[icut+1])
            rate_cut = rate[mass_mask]
            median = np.median(rate_cut)
            # Plot the histogram
            hist = axes[icut].hist(rate_cut, label=f'mass cut {mass_cuts[icut]}')[0]
            axes[icut].plot([median, median], [0, max(hist)], linestyle='dashed', color='red', 
                            label=f'median = {np.round(median, 3)}')
            axes[icut].legend(loc='best')
            
        plt.savefig(load_dir+f'snap_{snap_list[isnap][5:]}_cuts')
        plt.close()