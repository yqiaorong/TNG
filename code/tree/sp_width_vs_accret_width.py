import numpy as np
import os
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')

type = 'Hydro' # [DM / Hydro]
sim = 'MTNG' # [TNG300 / MTNG]
acc_width_type = 'percentile' # [std / percentile]


# Load accretion rate width
root_dir = f'result/accretion_rate_plot/{sim}/'
if sim == 'TNG300': 
    acc_data = np.load(f'{root_dir}/sim_205_1250_{type}/TNG300_{type}_accret_stats.npy', allow_pickle=True).item()
    snap_list = [8, 13, 25, 40, 67, 99]
elif sim == 'MTNG':
    acc_data = np.load(f'{root_dir}/{type}-Arepo/MTNG_{type}_accret_stats.npy', allow_pickle=True).item()
    snap_list = [264, 214, 151]
    
    
    
# Load all data
redshifts = acc_data['redshifts']
mass_cuts = acc_data['mass_cuts'] 
median_array = acc_data['accret_med']

if acc_width_type == 'std':
    width_array = acc_data['accret_std']
elif acc_width_type == 'percentile':
    lowp_array  = acc_data['accret_low']
    highp_array = acc_data['accret_high']
    width_array = highp_array - lowp_array
    del lowp_array, highp_array



# Select mass cuts with splashback features
if sim == 'MTNG':
    redshifts = redshifts[:len(snap_list)]
    
    cut_start_idx = np.where(mass_cuts == 3)[0][0]
    mass_cuts = mass_cuts[cut_start_idx:]
    
    median_array = median_array[:len(snap_list), cut_start_idx:]
    if acc_width_type == 'std':
        width_array = width_array[:len(snap_list), cut_start_idx:]
    elif acc_width_type == 'percentile':
        width_array = width_array[:len(snap_list), cut_start_idx:]
print('redshifts: ', np.round(redshifts, 3)) # From low z (present) to high z
print('mass cuts: ', mass_cuts)
print('(redshifts, mass_cuts) ', width_array.shape) # (redshifts, mass_cuts) 



# Load splashback features width
Nboots = 1024
feat_idx = 2 # width
boots_dir = f'result/bootstrap_stats/{sim}/'
if sim == 'TNG300': 
    sp_dir = f'{boots_dir}/sim_205_1250_{type}/Nboots_{Nboots}/'
elif sim == 'MTNG':
    sp_dir = f'{boots_dir}/{type}-Arepo/MTNG-L500-4320-A/output/Nboots_{Nboots}/'
TNG_list = os.listdir(sp_dir)

TNG_all_data, TNG_z = [], []
for snap in snap_list: # from low z (present) to high z 
    # Load data
    data = np.load(sp_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    
    z = data['z']
    TNG_z.append(z)
    
    data = data['final_results']
    TNG_all_data.append(data)
print('Sanity check of redshifts: ', np.round(TNG_z, 3))



# Set up the plot
sns.set(style="whitegrid")
fig, axs = plt.subplots(1, 1, dpi=500, figsize=(5, 4), sharey=True)
if type == 'Hydro':
    cmap_name = 'autumn'
elif type == 'DM':
    cmap_name = 'winter'
cmap = plt.get_cmap(cmap_name, len(TNG_list))


tot_plot_x, tot_plot_y = [], []
print(f'Number of mass cut {len(mass_cuts)-1}')
for cut_idx in range(len(mass_cuts)-1): # from low cut to high cut
    mask = median_array[:, cut_idx] != 0
    plot_x = width_array[mask, cut_idx]
    all_data = [d for d, m in zip(TNG_all_data, mask) if m==True]
    print('Number of snaps', len(all_data), f'in mass cut {mass_cuts[cut_idx]}')
    
    # Sort all data
    plot_y, plot_y_min, plot_y_max = [], [], []
    
    for snap_idx in range(len(all_data)): 
        print('(mass_cuts, feats, stats,)', all_data[snap_idx].shape, cut_idx)
        if all_data[snap_idx].shape[0] > cut_idx:
            plot_y.append(all_data[snap_idx][cut_idx, feat_idx, 1])
            plot_y_min.append(all_data[snap_idx][cut_idx, feat_idx, 0])
            plot_y_max.append(all_data[snap_idx][cut_idx, feat_idx, 2])
        else:
            pass
    
    axs.scatter(plot_x[np.argsort(plot_x)], np.array(plot_y)[np.argsort(plot_x)], 
                color=cmap(cut_idx / len(mass_cuts)),
                label=r'$10^{%.1f}$'%mass_cuts[cut_idx]+'~'
                +r'$10^{%.1f}$ '%mass_cuts[cut_idx+1]+'$M_\\odot$/h', 
                s=10)
    
    tot_plot_x.append(plot_x[np.argsort(plot_x)])
    tot_plot_y.append(np.array(plot_y)[np.argsort(plot_x)])
tot_plot_x = [item for sublist in tot_plot_x for item in sublist]
tot_plot_y = [item for sublist in tot_plot_y for item in sublist]



r_value, p_value = pearsonr(tot_plot_x, tot_plot_y)
if type == 'DM':
    reg_color = "blue"
elif type == 'Hydro':
    reg_color = "red"
axs = sns.regplot(x=tot_plot_x, y=tot_plot_y, ci=95, 
                  scatter=False, line_kws={"color": reg_color})
plt.title(f'r = {r_value:.2f}, p = {p_value:.3f}', loc='center')

# Final edit
axs.set_xlabel(f'Accretion rate width ({acc_width_type})')
axs.set_ylabel(f"Splashback feature (width)")
axs.legend()

# Save the plot
plt.savefig(f'{root_dir}/{sim}_{type}_sp_width_accret_width_{acc_width_type}')
plt.close()