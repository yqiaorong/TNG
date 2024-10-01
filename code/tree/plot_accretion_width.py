import os
import argparse
import numpy as np
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--z_or_Omz', default=None, type=str)
args = parser.parse_args()

# Set up the plot
fig, ax = plt.subplots(2, 2, figsize=(10, 6), sharex=True, sharey=True, 
                       # constrained_layout=False
                       )
# fig.subplots_adjust(hspace=0, wspace=0)
# plt.tight_layout(pad=1.0, h_pad=0, w_pad=0)

ax[0][0].set_ylabel('accretion rate width')
ax[1][0].set_ylabel('accretion rate width')
ax[1][1].set_xlabel(args.z_or_Omz)
ax[1][0].set_xlabel(args.z_or_Omz)        
ax[0][0].set_title(f'MTNG-Hydro')
ax[0][1].set_title(f'MTNG-DM')



# Run the scripts
# os.system('python3 code/tree/plot_accretion_hist2.py --sim DM')
# os.system('python3 code/tree/plot_accretion_hist2.py --sim Hydro')



def Om_z(z, Omega0, OmegaLambda):
    return Omega0*(1+z)**3/(Omega0*(1+z)**3 + OmegaLambda)

root_dir = f'result/accretion_rate_plot_new/'

#######################################################################################################
# MTNG-Hydro
#######################################################################################################

data = np.load(f'{root_dir}/Hydro-Arepo/MTNG-L500-4320-A/MTNG_Hydro_accret_stats.npy', 
               allow_pickle=True).item()

# Get Omega0(z)
Omega0 = data['Omega0']
print('Om:', Omega0)
OmegaLambda = data['OmegaLambda']
print('OLambda', OmegaLambda)
redshifts = data['redshifts']
Omega0z = Om_z(redshifts, Omega0, OmegaLambda)

redshifts = data['redshifts']
mass_cuts = data['mass_cuts']
std_array = data['accret_std']
median_array = data['accret_med']
lowp_array   = data['accret_low']
highp_array  = data['accret_high']

cmap = plt.get_cmap('autumn', len(mass_cuts))
for icut in range(len(mass_cuts)-1):
     mask = median_array[:, icut] != 0 # Remove zero terms

     if args.z_or_Omz == 'z':
          ax[0][0].plot(redshifts[mask], std_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                    label='std '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                         +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
          ax[1][0].plot(redshifts[mask], highp_array[mask, icut]-lowp_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                    label='width '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                         +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h', linestyle='--')  
     elif args.z_or_Omz == 'Omz':
          ax[0][0].plot(Omega0z[mask], std_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                    label='std '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                         +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
          ax[1][0].plot(Omega0z[mask], highp_array[mask, icut]-lowp_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                    label='width '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                         +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h', linestyle='--') 

#######################################################################################################
# MTNG-DM
#######################################################################################################

data = np.load(f'{root_dir}/DM-Arepo/MTNG-L500-4320-A/MTNG_DM_accret_stats.npy', 
               allow_pickle=True).item()

# Get Omega0(z)
Omega0 = data['Omega0']
print('Om:', Omega0)
OmegaLambda = data['OmegaLambda']
print('OLambda', OmegaLambda)
redshifts = data['redshifts']
Om_z_list = Om_z(redshifts, Omega0, OmegaLambda)

redshifts = data['redshifts']
mass_cuts = data['mass_cuts']
std_array = data['accret_std']
median_array = data['accret_med']
lowp_array   = data['accret_low']
highp_array  = data['accret_high']

cmap = plt.get_cmap('winter', len(mass_cuts))
for icut in range(len(mass_cuts)-1):
     mask = median_array[:, icut] != 0 # Remove zero terms
    
     if args.z_or_Omz == 'z':
          ax[0][1].plot(redshifts[mask], std_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                    label='std '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                         +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
          ax[1][1].plot(redshifts[mask], highp_array[mask, icut]-lowp_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                    label='width '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                         +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h', linestyle='--')
     elif args.z_or_Omz == 'Omz':
          ax[0][1].plot(Omega0z[mask], std_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                    label='std '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                         +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h')
          ax[1][1].plot(Omega0z[mask], highp_array[mask, icut]-lowp_array[mask, icut], color=cmap(icut/len(mass_cuts)), 
                    label='width '+r'$10^{%.1f}$'%mass_cuts[icut]+'~'
                         +r'$10^{%.1f}$ '%mass_cuts[icut+1]+'$M_\\odot$/h', linestyle='--')

# Save figure
for i in range(2):
    for j in range(2):
        ax[i][j].legend(loc='best')
plt.savefig(f'{root_dir}/MTNG_accretion_width_vs_{args.z_or_Omz}')
plt.close()