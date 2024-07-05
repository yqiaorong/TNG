from matplotlib import pyplot as plt 
plt.style.use('code/plot/style.mplstyle')
import h5py
import illustris_python as il
import numpy as np
import os
import math

boxsize = 205
res = 1250
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/' + 'L%dn%dTNG/output'%(boxsize,res)



snaps = [33, 40, 50, 67, 78, 99]
mass_cut = ['2', '2-5', '3', '3-5']
mass_cut_idx = 0



# Setup the plot
fig = plt.figure(figsize=(10,5), dpi=400)
gs = fig.add_gridspec(2, len(snaps), hspace=0, wspace=0)
axs = gs.subplots(sharex='col', sharey='row')


data_path = f'result/bootstrap'
save_dir = f'{data_path}/sim_{boxsize}_{res}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
    

# Redshift data
# z = []
for i, snap in enumerate(snaps):
    # Load redshifts
    with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        h = header['HubbleParam']
        
    # Load data
    data = np.load(data_path+f'/snap_{snap}_Rsp_stats.npy', allow_pickle = True).item()
    median_idx = data['median_idx_in_boots']
    
    # select median data
    profile = np.load(f'result/bootstrap/snap_{snap}/data'+
            f'/mass_cut_{mass_cut[mass_cut_idx]}/boots_{median_idx[mass_cut_idx]}.npy',
            allow_pickle=True).item()
    Rsp = data['final_results'][mass_cut_idx, 0, 1] * h / scale_factor / profile['R200_median'] # dimensionless
    del data

    # density profile
    ax0 = axs[0, i]
    ax0.errorbar(profile['radius'], profile['rho'], yerr = profile['rho_err'].T, color='r',
                 label=f'z={np.round(z, 3)}') # dimensionless radius
    ax0.axvline(x=Rsp, color='b', linestyle='--', linewidth=1)                              # dimensionless radius
    # gradient profile
    ax1 = axs[1, i]
    ax1.errorbar(profile['radius'], profile['slope'], yerr = profile['slope_err'].T, color='r',
                 label=f'z={np.round(z, 3)}') # dimensionless radius
    ax1.plot(profile['fitted_radius'], profile['fitted_slope'])                                 # dimensionless radius
    ax1.axvline(x=Rsp, color='b', linestyle='--', linewidth=1)                                  # dimensionless radius
    
    # general settings
    ax0.set_yscale('log')
    ax0.set_ylabel(r"Mass density$\rho$/$\rho_c$")
    ax0.legend(loc='best')
    ax1.set_ylabel(r'd log \rho / d \og r')
    ax1.set_ylim(-4,-0)
    ax1.legend(loc='best')

for ax in fig.get_axes():
    ax.set_xlabel(r"Radius r/$R_{200}$")
    ax.set_xscale('log')
    ax.label_outer()


plt.savefig(f'result/bootstrap/sim_{boxsize}_{res}/mass_cut_{mass_cut[mass_cut_idx]}_profiles')
plt.close
    
    

    
    
