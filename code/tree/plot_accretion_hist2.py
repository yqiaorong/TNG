import pandas as pd
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import h5py
import illustris_python as il
import numpy as np
import argparse
import os
from tqdm import tqdm

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim', default='DM', type=str)
args = parser.parse_args()

print('')
print(f'>>> Accretion rate histograms <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

simpath = f'{args.sim}-Arepo/MTNG-L500-4320-A/'
basePath = f'/virgotng/mpa/MTNG/{simpath}/'

load_dir = f'result/DMhalo_mass_table_new/{simpath}/'
save_dir = f'result/accretion_rate_plot_new/{simpath}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
snaps = [264, 214, 151, 94, 69, 51]

# Load df
h = 0.6774
mass_fname = "snap_{}_FPGrMass.npy"

accretion_df = np.load(load_dir+'accretion_rates.npy')



# Load redshift
redshifts, scale_factors = [], []
for snap in tqdm(snaps[:-1], desc='load z'):
    with h5py.File(il.snapshot.snapPath(basePath+'output', snap), 'r') as f:
        header = dict(f['Header'].attrs.items())

        Omega0 = 0.3089
        OmegaLambda = 0.6911
        
        a = header['Time']
        z = header['Redshift']
        
        redshifts.append(z)
        scale_factors.append(a)
redshifts     = np.array(redshifts)
scale_factors = np.array(scale_factors)



# 1. plot the total accretion rate per snapshots
tot_median, tot_std, tot_low_bound, tot_high_bound = [], [], [], []
for isnap, snap in enumerate(tqdm(snaps[:-1], desc='calc distri over snaps')):
    
    rate = accretion_df[:, isnap]
    rate = rate[rate != 0]
    # Sanity check
    if np.all(np.isfinite(rate) & (rate != 0)):
        print(f'No nan or inf in snap {snap}')
    
    median = np.median(rate)
    std    = np.std(rate)
    low_bound = np.percentile(rate, 16)
    high_bound = np.percentile(rate, 84)
    
    tot_median.append(median)
    tot_std.append(std)
    tot_low_bound.append(low_bound)
    tot_high_bound.append(high_bound)

    plt.figure()
    hist = plt.hist(rate, label=f'z = {np.round(redshifts[isnap], 2)}')[0]
    plt.plot([median, median], [0, max(hist)], linestyle='dashed', color='red', 
                label=f'median = {np.round(median, 3)}')
    plt.legend(loc='best')
    plt.savefig(save_dir+f'snap_{snap}')
        
tot_median     = np.array(tot_median)
tot_tsd        = np.array(tot_std)
tot_low_bound  = np.array(tot_low_bound)
tot_high_bound = np.array(tot_high_bound)



# 2. Plot the total accretion rate accross snapshots
plt.figure()
plt.plot(redshifts, tot_median, color='b', label='median')
plt.fill_between(redshifts, tot_median+tot_std, tot_median-tot_std, color='r', alpha=0.2, label='std')
plt.fill_between(redshifts, tot_low_bound, tot_high_bound, color='b', alpha=0.2, label='percentile')
plt.xlabel('z')
plt.ylabel('accretion rate')
plt.legend(loc='best')
plt.title(f'MTNG-{args.sim}')
plt.savefig(save_dir+f'tot_accretion_rate_vs_z_MTNG_{args.sim}')
plt.close()



# 3. Plot accretion rate per mass cut per snapshots
bin_start, bin_end = 1, 5
mass_cuts = np.arange(bin_start, bin_end, 0.5)
num_cuts  = int((bin_end-bin_start)/0.5)

median_array = np.zeros((len(snaps)-1, len(mass_cuts)-1))
std_array    = np.zeros((len(snaps)-1, len(mass_cuts)-1))
lowp_array   = np.zeros((len(snaps)-1, len(mass_cuts)-1))
highp_array  = np.zeros((len(snaps)-1, len(mass_cuts)-1))

for isnap in range(len(snaps)-1):
    
    fig, axes = plt.subplots(1, num_cuts, figsize=(15, 5))
    # Load mass
    mass = np.load(load_dir+mass_fname.format(snaps[isnap]))
    mass = mass / h
    
    # Exclude mass with no accretion rate
    rate = accretion_df[:, isnap]
    mask = rate != 0
    mass, rate = mass[mask], rate[mask]

    for icut, cut in enumerate(mass_cuts[:-1]): 
    # Assign the accretion rate between mass_cut m1 and m2 to the former m1
        # Get the mass bin
        cut_mask = (mass >= 10**cut) & (mass < 10**(cut+0.5))
        rate_cut = rate[cut_mask]
        
        if rate_cut.size == 0:
            pass
        else:
            # At later snaps, two normal distri of accret rates appear per mass cut per snap,
            # Can check the histogram per mass cut per snap. 
            # Here we just simply remove the higher accret rates distri.
            
            # rate_mask = (rate_cut <= 15)
            # rate_cut  = rate_cut[rate_mask]
            
            if len(rate_cut) < 30:
                pass
            else:
                median     = np.median(rate_cut)
                std        = np.std(rate_cut)
                low_bound  = np.percentile(rate_cut, 16)
                high_bound = np.percentile(rate_cut, 84)
                
                # Plot the histogram
                hist = axes[icut].hist(rate_cut, label=f'mass cut {cut} ~ {cut+0.5}')[0]
                axes[icut].plot([median, median], [0, max(hist)], linestyle='dashed', color='red', 
                                label=f'median = {np.round(median, 3)}')
                axes[icut].legend(loc='best')
                
                median_array[isnap, icut] = median
                std_array[isnap, icut]    = std
                lowp_array[isnap, icut]   = low_bound
                highp_array[isnap, icut]  = high_bound
        
    plt.savefig(save_dir+f'snap_{snaps[isnap]}_cuts')
    plt.close()
        


# 4. Plot accretion rate per mass cut across snapshots
fig, ax = plt.subplots(1, 1)
if args.sim == 'DM':
    cmap_name = 'winter'
elif args.sim == 'Hydro':
    cmap_name = 'autumn'
cmap = plt.get_cmap(cmap_name, len(mass_cuts))

for icut in range(len(mass_cuts)-1):
    mask = median_array[:, icut] != 0 # Remove zero terms
    ax.plot(redshifts[mask], median_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
             label=r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                  +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
    # plt.fill_between(redshifts[mask], median_array[mask, icut]+std_array[mask, icut],
    #                             median_array[mask, icut]-std_array[mask, icut], alpha=0.2)
    plt.fill_between(redshifts[mask], lowp_array[mask, icut], highp_array[mask, icut], 
                     alpha=0.1, color=cmap(icut/len(mass_cuts)))
ax.set_xlabel('z')
ax.set_ylabel('accretion rate')
ax.legend(loc='best')
ax.set_title(f'MTNG-{args.sim}')
plt.savefig(save_dir+f'accretion_rate_vs_z_MTNG_{args.sim}')
plt.close()

# Save accret per mass cut per snap
save_dict = {'scale_factors':scale_factors, 'redshifts': redshifts, 'Omega0': Omega0, 'OmegaLambda': OmegaLambda,
             'mass_cuts': mass_cuts,
             'accret_std': std_array,  
             'accret_med': median_array,
             'accret_low': lowp_array, 
             'accret_high': highp_array}
np.save(f'{save_dir}/MTNG_{args.sim}_accret_stats', save_dict)