import numpy as np
import os

sim = 'Hydro'
load_dir = f'result/DMhalo_density_profiles_phys/MTNG/{sim}-Arepo/MTNG-L500-4320-A/'
# load_dir = f'result/DMhalo_density_profiles_phys/TNG300/sim_205_1250_{sim}/'
print(sim, load_dir)

for snap in os.listdir(load_dir):
    file_dir = load_dir+snap+'/final_densities/'
    for file in os.listdir(file_dir):
        f = np.load(file_dir+file, allow_pickle=True).item()
        z = f['z']
        print(snap, np.round(z,3))
        halo_M_Mean200 = f['halo_M_Mean200']
        bin_start, bin_end = int(file[:-4].split('-')[1]), int(file[:-4].split('-')[2])
        # if bin_end - bin_start != 5:
        #     mass_bins = np.linspace(bin_start, bin_end, num=int((bin_end-bin_start)/5)+1)
        #     # count number of elements between mass bins
        #     for b in mass_bins[:-1]:
        #         print(snap, b, b+5, np.where((10**(b/10) <= halo_M_Mean200) & 
        #                                 (halo_M_Mean200 < 10**((b+5)/10)))[0].shape
        #               )
        # else:
        #     print(snap, file, halo_M_Mean200.shape)