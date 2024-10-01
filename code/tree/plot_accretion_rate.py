import os
import numpy as np
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')

fig, ax = plt.subplots(1, 2, figsize=(8, 4))

root_dir = f'result/accretion_rate_plot_new/'

# Run the scripts
# os.system('python3 code/tree/plot_accretion_hist2.py --sim DM')
# os.system('python3 code/tree/plot_accretion_hist2.py --sim Hydro')

# MTNG-Hydro
data = np.load(f'{root_dir}/Hydro-Arepo/MTNG-L500-4320-A/MTNG_Hydro_accret_stats.npy', 
               allow_pickle=True).item()

redshifts = data['redshifts']
mass_cuts = data['mass_cuts']
std_array = data['accret_std']
median_array = data['accret_med']
lowp_array   = data['accret_low']
highp_array  = data['accret_high']

cmap = plt.get_cmap('autumn', len(mass_cuts))
for icut in range(len(mass_cuts)-1):
    
    mask = median_array[:, icut] != 0 # Remove zero terms
    ax[0].plot(redshifts[mask], median_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                label=r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                    +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
    ax[0].fill_between(redshifts[mask], lowp_array[mask, icut], highp_array[mask, icut], 
                       alpha=0.1, color=cmap(icut/len(mass_cuts)))

#TNG300_DM
data = np.load(f'{root_dir}/DM-Arepo/MTNG-L500-4320-A/MTNG_DM_accret_stats.npy', 
               allow_pickle=True).item()

redshifts = data['redshifts']
mass_cuts = data['mass_cuts']
std_array = data['accret_std']
median_array = data['accret_med']
lowp_array   = data['accret_low']
highp_array  = data['accret_high']

cmap = plt.get_cmap('winter', len(mass_cuts))
for icut in range(len(mass_cuts)-1):
    
    mask = median_array[:, icut] != 0 # Remove zero terms
    ax[1].plot(redshifts[mask], median_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                label=r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                    +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
    ax[1].fill_between(redshifts[mask], lowp_array[mask, icut], highp_array[mask, icut], 
                     alpha=0.1, color=cmap(icut/len(mass_cuts)))
    
# General settings and save fig
ax[0].set_ylabel('accretion rate width')
ax[0].set_xlabel('z')
ax[1].set_xlabel('z')
ax[0].legend(loc='best')
ax[1].legend(loc='best')
ax[0].set_title(f'MTNG-Hydro')
ax[1].set_title(f'MTNG-DM')
plt.savefig(f'{root_dir}/MTNG_accretion_rate_vs_z')
plt.close()