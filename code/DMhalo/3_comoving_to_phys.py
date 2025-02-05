"""WARNINGS: Commmented parts are finished which should never to run again.
   All the data computed in this script is saved in result/DMhalo_density_profiles_phys"""

import os
import numpy as np
from astropy.cosmology import Planck15



###########################################################################################
# Convwet the comoving data to physical data
###########################################################################################

# sim = 'MTNG/Hydro-Arepo/MTNG-L500-4320-A/'
# snaps = [129, 151, 179, 214, 237, 264]

# sim = 'TNG300/sim_205_1250_Hydro/'
# snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]

# for snap in snaps:
    
#     load_dir = f'result/DMhalo_density_profiles_raw2/{sim}/snap_{snap}/final_densities/'
#     save_dir = f'result/DMhalo_density_profiles_phys/{sim}/snap_{snap}/final_densities/'
#     if not os.path.exists(save_dir):
#         os.makedirs(save_dir)
    
#     # Load the data
#     for fname in os.listdir(load_dir):
#         print(fname)
#         data = np.load(load_dir + fname, allow_pickle=True).item()
#         print(data.keys())
        
#         h, scale_factor, z = data['h'], data['scale_factor'], data['z']
        
#         # Change mass
#         halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 Msun/h]
#         halo_M_Mean200p = halo_M_Mean200 / h    # [10^10 Msun]
#         del halo_M_Mean200
        
#         # Change radius
#         halo_R_Mean200 = data['halo_R_Mean200']             # [ckpc/h]
#         halo_R_Mean200p = halo_R_Mean200 * scale_factor / h # [kpc]
#         del halo_R_Mean200
        
#         radial_bins = data['radial_bins']             # [ckpc/h]
#         radial_binsp = radial_bins * scale_factor / h # [kpc]
#         del radial_bins
        
#         # Change density
#         densities = data['densities']                       # [(Msun/h) / (ckpc/h)^3]
#         densitiesp = densities * (h**2) / (scale_factor**3) # [(Msun) / (kpc)^3]
#         del densities
        
#         del data
        
#         # Calculate the critical density
#         rho_c = Planck15.critical_density(z).to('Msun/kiloparsec**3') # [(Msun) / (kpc)^3]
#         print(rho_c)
        
#         # Save the new physical quantities
#         save_dict = {'halo_R_Mean200': halo_R_Mean200p, # [kpc]
#                     'halo_M_Mean200': halo_M_Mean200p, # [10^10 Msun]
#                     'h': h, 'scale_factor': scale_factor, 'z': z,
#                     'radial_bins': radial_binsp,       # [kpc]
#                     'densities': densitiesp,           # [(Msun) / (kpc)^3]
#                     'rho_c': rho_c                     # [(Msun) / (kpc)^3]
#                     }    
#         np.save(save_dir + fname, save_dict)
#         print(fname + ' saved')
#         print('')
        
        
        
###########################################################################################
# Convwet MTNG data from Mpc to kpc
###########################################################################################

# sim = 'MTNG/Hydro-Arepo/MTNG-L500-4320-A/'
# sim = 'MTNG/DM-Arepo/MTNG-L500-4320-A/'
# snaps = [129, 151, 179, 214, 237, 264]

# for snap in snaps:
    
#     load_dir = f'result/DMhalo_density_profiles_phys/{sim}/snap_{snap}/final_densities/'
    
#     # Load the data
#     for fname in os.listdir(load_dir):
#         print(fname)
#         data = np.load(load_dir + fname, allow_pickle=True).item()
#         print(data.keys())
        
        # data['halo_R_Mean200'] = data['halo_R_Mean200']*1000 # [kpc]
        # data['radial_bins'] = data['radial_bins']*1000       # [kpc]
        # data['densities'] = data['densities'] / (10**9)      # [(Msun) / (kpc)^3]
        
        # # Save the new physical quantities
        # np.save(load_dir + fname, data)

        # print(fname + ' saved')
        # print('')