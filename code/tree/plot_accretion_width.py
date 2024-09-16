import os
import numpy as np
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')

fig, ax = plt.subplots(1, 2, figsize=(10, 4))

root_dir = 'result/accretion_rate_plot/'

# Run the scripts
# os.system('python3 code/tree/plot_accretion_hist.py --sim DM')
# os.system('python3 code/tree/plot_accretion_hist.py --sim Hydro')

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
     
     print(mass_cuts[icut])
     mask = median_array[:, icut] != 0 # Remove zero terms
     print(std_array[mask, icut])
     print(highp_array[mask, icut]-lowp_array[mask, icut])
     print('')
     ax[0].plot(redshifts[mask], std_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
          label='std '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
               +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
     ax[0].plot(redshifts[mask], highp_array[mask, icut]-lowp_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
          label='width '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
               +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h', linestyle='--')  

# MTNG-DM
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
    
    ax[1].plot(redshifts[mask], std_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
             label='std '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                  +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
    ax[1].plot(redshifts[mask], highp_array[mask, icut]-lowp_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
             label='width '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                  +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h', linestyle='--')
    
# General settings and save fig
ax[0].set_ylabel('accretion rate width')
ax[0].set_xlabel('z')
ax[1].set_xlabel('z')
ax[0].legend(loc='best')
ax[1].legend(loc='best')
ax[0].set_title(f'MTNG-Hydro')
ax[1].set_title(f'MTNG-DM')
plt.savefig(f'{root_dir}/accretion_width_vs_z')
plt.close()