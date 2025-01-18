"""The density profiles."""
import matplotlib.pyplot as plt
plt.style.use('code/style.mplstyle')
import numpy as np
import os

# =============================================================================
# Load sample data
# =============================================================================

root_dir = 'result/bootstrap_stats_phys/'
sim_dir = f'TNG300/sim_205_1250_Hydro/'

# Load stat data
data = np.load(f'{root_dir}/{sim_dir}/Nboots_1024/snap_99_Rsp_stats.npy', allow_pickle = True).item()

z = data['z']
h = data['h']

mass_bins = data['mass_bins']
# Select halos with mass 10^13 MSun
cut = 2
mass_cut_idx = np.where(mass_bins == cut)[0][0]

median_idx = data['median_idx_in_boots']
phy_Rsp = data['final_results'][mass_cut_idx, 0, 1] # [kpc]

# Load profile data
profile_dir = 'result/bootstrap_phys/'
profile = np.load(f'{profile_dir}/{sim_dir}/snap_99/Nboots_1024/data/'+
                    f'/mass_cut_{cut}/boots_{median_idx[mass_cut_idx]}.npy',
                    allow_pickle=True).item()

phy_R200 = profile['R200_median']                   # [kpc]
scale_Rsp =  phy_Rsp / phy_R200                     # [dimensionless]
del data

# Compute softening length in R200
# The gravitational softening length (after multiplied with 2.8):
# TNG300: 4kpc
scale_Rsoft = 4 / phy_R200              # [dimensionless]

# =============================================================================
# Plot
# =============================================================================
   
fig, axs = plt.subplots(2, 1, figsize=(4,6), dpi=500, sharex=True, constrained_layout=True)

# Color
from matplotlib.colors import LinearSegmentedColormap
red_list = ['#EE9D9F', '#DE6A69', '#C84747', '#982B2D','#6A0624', '#3D011A'] # light to dark
blue_list = ['#89CAEA','#4596CD','#0B75B3','#015696','#012A61','#053061'] # light to dark
green_list = ['#C8D7B4', '#A5BD8A','#87A368','#739353','#50673A','#3D4F2F','#333F29','#2C3725'] # light to dark

# density profile
axs[0].errorbar(profile['radius'], profile['rho'], yerr = profile['rho_err'].T, 
                color='#2C3725', fmt='.', markersize=3)
axs[0].plot(profile['fitted_radius'], profile['fitted_rho'], color='#50673A')  # dimensionless radius
axs[0].axvline(x=scale_Rsp, color='#50673A', linestyle='dashed', linewidth=1)                              # dimensionless radius
axs[0].axvline(x=scale_Rsoft, color='#50673A', linestyle='dotted', linewidth=1)                              # dimensionless radius

axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_ylabel(r'$\rho / \rho_c$')

# gradient profile
axs[1].errorbar(profile['radius'], profile['slope'], yerr = profile['slope_err'].T,
                color='#2C3725', fmt='.', markersize=3)
axs[1].plot(profile['fitted_radius'], profile['fitted_slope'], color='#50673A')  # dimensionless radius
axs[1].axvline(x=scale_Rsp, color='#50673A', linestyle='dashed', linewidth=1)                              # dimensionless radius
axs[1].axvline(x=scale_Rsoft, color='#50673A', linestyle='dotted', linewidth=1)                              # dimensionless radius

axs[1].set_xscale('log')
axs[1].set_ylabel(r'$d \log \rho / d \log r$')
axs[1].set_xlabel(r'Radius $r / R_{200}$')

# =============================================================================
# Annotate depth and width
# =============================================================================

# Annotate depth (vertical arrow)
min_grad = np.min(profile['fitted_slope'])
min_grad_idx = np.argmin(profile['fitted_slope'])
min_grad_R = profile['fitted_radius'][min_grad_idx]

left_data = profile['fitted_slope'][:min_grad_idx]
right_data = profile['fitted_slope'][min_grad_idx:]
max_grad = np.max(right_data)

axs[1].annotate('', xy=(min_grad_R, min_grad), xytext=(min_grad_R, max_grad),
                arrowprops=dict(arrowstyle='<->', color=red_list[-3], lw=1.5))
# axs[1].text(1.05, -2, 'Depth', color='green', fontsize=10)

# Annotate width (horizontal arrow)
depth = max_grad - min_grad
half_grad = min_grad + depth/2

left_idx = np.argmin(np.abs(left_data - half_grad))
right_idx = min_grad_idx + np.argmin(np.abs(right_data - half_grad))

left_R = profile['fitted_radius'][left_idx]
right_R = profile['fitted_radius'][right_idx]

axs[1].annotate('', xy=(left_R, half_grad), xytext=(right_R, half_grad),
                arrowprops=dict(arrowstyle='<->', color=blue_list[-3], lw=1.5))
# axs[1].text(1.1, -2.5, 'Width', color='blue', fontsize=10)

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig1.png')
plt.close()