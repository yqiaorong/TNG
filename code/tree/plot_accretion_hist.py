import pandas as pd
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import h5py
import illustris_python as il
import numpy as np
import argparse
import os

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM', default='', type=str)
args = parser.parse_args()

print('')
print(f'>>> Plot accretion rate histogram <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

DM = args.DM
boxsize, res = 205, 1250

data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if args.DM == 'DM':
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'_DM/output/'
elif args.DM == 'Hydro':
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'/output/'

load_dir = f'result/DMhalo_mass_table/sim_{boxsize}_{res}_{DM}/'
save_dir = f'result/accretion_rate_plot/TNG300/sim_{boxsize}_{res}_{DM}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
    
    
# Load df
mass_df      = pd.read_csv(load_dir+'mass_table.csv', index_col=0)
accretion_df = pd.read_csv(load_dir+'accretion_table.csv', index_col=0)

# Select only accross cosmological time
mass_df      = mass_df.loc[['snap_8', 'snap_13', 'snap_25', 'snap_40', 'snap_67', 'snap_99']]
accretion_df = accretion_df.loc[['snap_8', 'snap_13', 'snap_25', 'snap_40', 'snap_67', 'snap_99']]

snap_list = mass_df.index
print(snap_list)


# Load redshift
redshifts, scale_factors = [], []
for snap in snap_list:
    with h5py.File(il.snapshot.snapPath(basePath, int(snap[5:])), 'r') as f:
        header = dict(f['Header'].attrs.items())
        
        Omega0 = header['Omega0']
        OmegaLambda = header['OmegaLambda']
        
        z = header['Redshift']
        a = header['Time']
        
        redshifts.append(z)
        scale_factors.append(a)
redshifts     = np.array(redshifts)
scale_factors = np.array(scale_factors)



# plot the total accretion rate per snapshots
tot_median, tot_std, tot_low_bound, tot_high_bound = [], [], [], []
for isnap in range(len(snap_list)):
    if isnap != 0:
        rate = accretion_df.iloc[isnap, :]
        rate = rate.replace([np.inf, -np.inf], np.nan)
        rate = rate.dropna().to_numpy()
        
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
        plt.savefig(save_dir+f'snap_{snap_list[isnap][5:]}')
        
tot_median     = np.array(tot_median)
tot_tsd        = np.array(tot_std)
tot_low_bound  = np.array(tot_low_bound)
tot_high_bound = np.array(tot_high_bound)



# Plot the total accretion rate accross snapshots
plt.figure()
plt.plot(redshifts[1:], tot_median, color='b', label='median')
plt.fill_between(redshifts[1:], tot_median+tot_std, tot_median-tot_std, color='r', alpha=0.2, label='std')
plt.fill_between(redshifts[1:], tot_low_bound, tot_high_bound, color='b', alpha=0.2, label='percentile')
plt.xlabel('z')
plt.ylabel('accretion rate')
plt.legend(loc='best')
plt.title(f'TNG300{DM}')
plt.savefig(save_dir+f'tot_accretion_rate_vs_z_TNG300_{DM}')
plt.close()



# Plot accretion rate per mass cut per snapshots
mass_cuts = np.arange(1, 4.5, 0.5)
num_cuts  = int((4.5-1)/0.5)

median_array = np.zeros((len(snap_list), len(mass_cuts)-1))
std_array    = np.zeros((len(snap_list), len(mass_cuts)-1))
lowp_array   = np.zeros((len(snap_list), len(mass_cuts)-1))
highp_array  = np.zeros((len(snap_list), len(mass_cuts)-1))

for isnap in range(len(snap_list)):
    
    if isnap != 0: # Assign the accretion rate between snapshots s1 and s2 to the later s2
        fig, axes = plt.subplots(1, num_cuts, figsize=(15, 5))
        
        mass = mass_df.iloc[isnap, :]
        rate = accretion_df.iloc[isnap, :]
        
        # Drop inf and nan
        rate = rate.replace([np.inf, -np.inf], np.nan)
        valid_mask = ~np.isnan(rate)
        rate, mass = rate[valid_mask], mass[valid_mask]

        for icut, cut in enumerate(mass_cuts[:-1]): 
        # Assign the accretion rate between mass_cut m1 and m2 to the former m1
            # Get the mass bin
            mass_mask = (mass >= 10**mass_cuts[icut]) & (mass < 10**mass_cuts[icut+1])
            rate_cut = rate[mass_mask]
            
            if rate_cut.empty or len(rate_cut) == 1:
                pass
            else:
                # At later snaps, two normal distri of accret rates appear per mass cut per snap,
                # Can check the histogram per mass cut per snap. 
                # Here we just simply remove the higher accret rates distri.
                rate_mask = (rate_cut <= 15)
                rate_cut  = rate_cut[rate_mask]
                
                median     = np.median(rate_cut)
                std        = np.std(rate_cut)
                low_bound  = np.percentile(rate_cut, 16)
                high_bound = np.percentile(rate_cut, 84)
                
                # Plot the histogram
                hist = axes[icut].hist(rate_cut, label=f'mass cut {mass_cuts[icut]}')[0]
                axes[icut].plot([median, median], [0, max(hist)], linestyle='dashed', color='red', 
                                label=f'median = {np.round(median, 3)}')
                axes[icut].legend(loc='best')
                
                median_array[isnap, icut] = median
                std_array[isnap, icut]    = std
                lowp_array[isnap, icut]   = low_bound
                highp_array[isnap, icut]  = high_bound
            
        plt.savefig(save_dir+f'snap_{snap_list[isnap][5:]}_cuts')
        plt.close()
        


# Plot accretion rate per mass cut across snapshots
fig, ax = plt.subplots(1, 1)
cmap = plt.get_cmap('autumn', len(mass_cuts))
for icut in range(len(mass_cuts)-1):
    
    mask = median_array[:, icut] != 0 # Remove zero terms
    ax.plot(redshifts[mask], median_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
             label=r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                  +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
    plt.fill_between(redshifts[mask], lowp_array[mask, icut], highp_array[mask, icut], 
                     alpha=0.1, color=cmap(icut/len(mass_cuts)))
ax.set_xlabel('z')
ax.set_ylabel('accretion rate')
ax.legend(loc='best')
ax.set_title(f'TNG300{DM}')
plt.savefig(save_dir+f'accretion_rate_vs_z_TNG300_{DM}')
plt.close()



# Save accret per mass cut per snap
save_dict = {'scale_factors':scale_factors, 'redshifts': redshifts, 'Omega0': Omega0, 'OmegaLambda': OmegaLambda,
             'mass_cuts': mass_cuts,
             'accret_std': std_array,  
             'accret_med': median_array,
             'accret_low': lowp_array, 
             'accret_high': highp_array}
np.save(f'{save_dir}/TNG300{DM}_accret_stats', save_dict)