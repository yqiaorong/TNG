from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import h5py
import illustris_python as il
import numpy as np
import os
import math

boxsize = 205
res = 1250
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/' + 'L%dn%dTNG/output'%(boxsize,res)

# Input
Nboots = 128


snaps = [33, 40, 50, 67, 78, 99]
mass_cut = ['1', '1.5', '2', '2-5', '3', '3-5']
mass_cut_idx = 2 # input



# Setup the plot
fig = plt.figure(figsize=(10,5), dpi=400)
gs = fig.add_gridspec(2, len(snaps), hspace=0, wspace=0)
axs = gs.subplots(sharex='col', sharey='row')

stats_path = f'result/bootstrap_stats/sim_{boxsize}_{res}/Nboots_{Nboots}/'  
    
### Calculate Msp

Msp_list, Rsp_list = [], []
for i, snap in enumerate(snaps):
    print(snap)
        
    # Load data
    data = np.load(stats_path+f'/snap_{snap}_Rsp_stats.npy', allow_pickle = True).item()
    z = data['z']
    scale_factor = 1 / (1+z)
    h = data['h']
    median_idx = data['median_idx_in_boots']

    # select median data
    profile = np.load(f'result/bootstrap/sim_205_1250/Nboots_{Nboots}/snap_{snap}/data'+
            f'/mass_cut_{mass_cut[mass_cut_idx]}/boots_{median_idx[mass_cut_idx]}.npy',
            allow_pickle=True).item()
    
    ### Rsp in [ckpc/h]
    from unyt import G, second, megaparsec, km, kiloparsec
    Rsp_phys = data['final_results'][mass_cut_idx, 0, 1]
    Rsp = Rsp_phys * h / scale_factor * kiloparsec               # [ckpc/h]
    
    ### Convert the radius from dimensionless to physical [ckpc/h]
    R200_median = profile['R200_median']                         # [ckpc/h]
    radius = profile['fitted_radius'] * R200_median * kiloparsec # [ckpc/h]
    del data
    
    ### Convert rho from dimensionless to [(MSun/h)/(ckpc/h)**3] ###
    
    # Compute the critical density
    Hubble = h * 100 * km / megaparsec / second
    rho_c = 3 * Hubble**2 / (8*np.pi*G) 
    rho_c.convert_to_units('Msun/kiloparsec**3') # [MSun/(kpc)**3]
    rho_c = rho_c * scale_factor**3 / h**2       # [(MSun/h)/(ckpc/h)**3]
    
    densities = profile['fitted_rho'] * rho_c    # [(MSun/h)/(ckpc/h)**3]

    ### Calculate Msp ###
    for irho, rho in enumerate(densities):
        r = radius[irho]
        if r > Rsp:
            pass
        else:
            if irho == 0:
                mass = 4/3 * np.pi * r**3 * rho
            else:
                dr = r - radius[irho-1] 
                delta_mass = 4 * np.pi * r**2 * dr * rho
                mass = mass + delta_mass
    print(f'halo mass Msp = {mass:.3f}/h and Rsp = {Rsp_phys:.3f} kpc at z = {z:.3f}')
    del mass