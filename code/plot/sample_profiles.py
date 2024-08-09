from matplotlib import pyplot as plt 
plt.style.use('code/plot/style.mplstyle')
import h5py
import illustris_python as il
import numpy as np
import os

boxsize = 205
res = 1250
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/' + 'L%dn%dTNG/output'%(boxsize,res)



snap = 33 # input
mass_cut = ['1', '1.5', '2', '2-5', '3', '3-5']
mass_cut_idx = 2 # input



# Setup the plot
fig = plt.figure(figsize=(2.5,4), 
                 dpi=400)
gs = fig.add_gridspec(2, 1, hspace=0, wspace=0)
axs = gs.subplots(sharex='col', sharey='row')


data_path = f'result/bootstrap/sim_{boxsize}_{res}'
save_dir = data_path
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
    

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
profile = np.load(f'{data_path}/snap_{snap}/data'+
        f'/mass_cut_{mass_cut[mass_cut_idx]}/boots_{median_idx[mass_cut_idx]}.npy',
        allow_pickle=True).item()
Rsp = data['final_results'][mass_cut_idx, 0, 1] * h / scale_factor / profile['R200_median'] # dimensionless
del data

# density profile
ax0 = axs[0]
ax0.errorbar(profile['radius'], profile['rho'], yerr = profile['rho_err'].T, color='r',
                label=f'z={np.round(z, 3)}') # dimensionless radius
ax0.axvline(x=Rsp, color='b', linestyle='--', linewidth=1)                              # dimensionless radius
# gradient profile
ax1 = axs[1]
ax1.errorbar(profile['radius'], profile['slope'], yerr = profile['slope_err'].T, color='r',
                label=f'z={np.round(z, 3)}') # dimensionless radius
ax1.plot(profile['fitted_radius'], profile['fitted_slope'])                             # dimensionless radius
ax1.axvline(x=Rsp, color='b', linestyle='--', linewidth=1)                              # dimensionless radius

# depth
grad_max = profile['fitted_slope'][-1]
grad_min = min(profile['fitted_slope'])
depth = grad_max - grad_min
radius_min_idx = np.argmin(profile['fitted_slope'])
radius_min = profile['fitted_radius'][radius_min_idx]

# ax1.annotate('', xy=(radius_min, grad_max), xytext=(radius_min, grad_min),
#              arrowprops=dict(arrowstyle='|-|', color='g', lw=1.5))

# width
half_grad = grad_min + depth/2
left_grads, right_grads = profile['fitted_slope'][:radius_min_idx], profile['fitted_slope'][radius_min_idx:]
left_idx = np.argmin(np.abs(left_grads - half_grad))
right_idx = radius_min_idx + np.argmin(np.abs(right_grads - half_grad))

# ax1.annotate('', xy=(profile['fitted_radius'][left_idx], profile['fitted_slope'][left_idx]), 
#                  xytext=(profile['fitted_radius'][right_idx], profile['fitted_slope'][right_idx]),
#              arrowprops=dict(arrowstyle='|-|', color='orange', lw=1.5))

# general settings
ax0.set_yscale('log')
ax0.set_ylabel(r"$\rho$/$\rho_c$")
ax0.legend(loc='best')
ax1.set_ylabel(r'd log $\rho$ / d \og r')
ax1.set_ylim(-4,-0)
ax1.legend(loc='best')

for ax in fig.get_axes():
    ax.set_xlabel(r"Radius r/$R_{200}$")
    ax.set_xscale('log')
    ax.label_outer()

plt.savefig(f'result/bootstrap/sim_{boxsize}_{res}/plot/'+
            f'sample_profile_snap{snap}_masscut{mass_cut[mass_cut_idx]}')
plt.close