save_root_dir = 'DMhalo_density_profiles'
sim = 'Hydro-Arepo/MTNG-L500-4320-A'
snapnum = 264
parent_dir = f'/nfs/mvogelsblab001/Users/s_qyu/{save_root_dir}/{sim}/snap_{snapnum}/'
save_dir = f'{parent_dir}/final_densities/'

import numpy as np
data = np.load(f'{save_dir}/bin-30-40.npy', allow_pickle=True).item() 
sample_idx = 10
sample_densities = data['densities'][sample_idx] / (10**9)
radial_bins = data['radial_bins'][sample_idx]

import matplotlib.pyplot as plt
plt.figure(figsize=(8, 6))
plt.plot(radial_bins, sample_densities, marker='o', linestyle='-')
plt.xscale('log')
plt.yscale('log')
plt.savefig('sample_profile')