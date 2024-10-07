import numpy as np
import os
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')

sim_type = 'Hydro' # [DM / Hydro]
sim = 'TNG300'       # [TNG300 / MTNG]
acc_width_type = 'percentile' # [std / percentile]


# Load accretion rate width
root_dir = f'result/accretion_rate_plot/{sim}/'
if sim == 'TNG300': 
    acc_data = np.load(f'{root_dir}/sim_205_1250_{sim_type}/TNG300_{sim_type}_accret_stats.npy', allow_pickle=True).item()
    snap_list = [8, 13, 25, 40, 67, 99]
elif sim == 'MTNG':
    acc_data = np.load(f'{root_dir}/{sim_type}-Arepo/MTNG_{sim_type}_accret_stats.npy', allow_pickle=True).item()
    snap_list = [151, 214, 264]
    
    
    
# Load all data
acc_z = np.round(acc_data['redshifts'], 3)  # From high z to low z (present)
acc_cut = acc_data['mass_cuts'][:-1] # from low cut to high cut
acc_med = acc_data['accret_med']     # (snaps, cuts,)

if acc_width_type == 'std':
    acc_width = acc_data['accret_std']   
elif acc_width_type == 'percentile':
    lowp_array  = acc_data['accret_low']
    highp_array = acc_data['accret_high']
    acc_width = highp_array - lowp_array
    del lowp_array, highp_array

# Select mass cuts with splashback features
# if sim == 'MTNG':
#     redshifts = redshifts[:len(snap_list)]
    
#     cut_start_idx = np.where(mass_cuts == 3)[0][0]
#     mass_cuts = mass_cuts[cut_start_idx:]
    
#     median_array = median_array[:len(snap_list), cut_start_idx:]
#     if acc_width_type == 'std':
#         width_array = width_array[:len(snap_list), cut_start_idx:]
#     elif acc_width_type == 'percentile':
#         width_array = width_array[:len(snap_list), cut_start_idx:]

print('redshifts: ', acc_z) # From low z (present) to high z
print('mass cuts: ', acc_cut)
print('(redshifts, mass_cuts) ', acc_width.shape) # (redshifts, mass_cuts) 



# Load splashback features width
Nboots = 1024
feat_idx = 2 # width
boots_dir = f'result/bootstrap_stats/{sim}/'
if sim == 'TNG300': 
    sp_dir = f'{boots_dir}/sim_205_1250_{sim_type}/Nboots_{Nboots}/'
elif sim == 'MTNG':
    sp_dir = f'{boots_dir}/{sim_type}-Arepo/MTNG-L500-4320-A/Nboots_{Nboots}/'
TNG_list = os.listdir(sp_dir)

TNG_data, TNG_z, TNG_cut = [], [], []
for snap in snap_list: # From low z (present) to high z
    # Load data
    data = np.load(sp_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    TNG_z.append(np.round(data['z'],3))
    TNG_cut.append(data['mass_bins'])
    TNG_data.append(data['final_results'])
print('Sanity check of redshifts: ', TNG_z)



# Set up the plot
sns.set(style="whitegrid")
fig, axs = plt.subplots(1, 1, dpi=500, figsize=(5, 4), sharey=True)
if sim_type == 'Hydro':
    cmap_name = 'autumn'
elif sim_type == 'DM':
    cmap_name = 'winter'
cmap = plt.get_cmap(cmap_name, len(TNG_list))


tot_plot_x, tot_plot_y = [], []
print(f'Number of mass cut {len(acc_cut)}')
for current_cut in acc_cut: # from low cut to high cut
    
    print('Current mass cut: ', current_cut)
    acc_cut_idx = np.where(acc_cut == current_cut)[0][0]
    acc_mask = acc_med[:, acc_cut_idx] != 0 # At specific mass cut, the snap values
    
    plot_x_width = acc_width[acc_mask, acc_cut_idx] # The width accross snaps per mass cut
    plot_x_z = acc_z[acc_mask]
    
    data_snap_idx = [int(np.where(TNG_z == z)[0][0]) for z in plot_x_z]
    select_data = [TNG_data[idx] for idx in data_snap_idx]
    select_cut  = [TNG_cut[idx] for idx in data_snap_idx]
    select_z    = [TNG_z[idx] for idx in data_snap_idx]
    
    
    # Sort all data
    plot_y, plot_y_min, plot_y_max, plot_y_z = [], [], [], []
    for z, cut, data in zip(select_z, select_cut, select_data): 
        print('(mass_cuts, feats, stats,)', data.shape, cut.shape, cut)
        try:
            data_cut_idx = np.where(cut == current_cut)[0][0]
            print(data_cut_idx)
            plot_y.append(data[data_cut_idx, feat_idx, 1])
            plot_y_min.append(data[data_cut_idx, feat_idx, 0])
            plot_y_max.append(data[data_cut_idx, feat_idx, 2])
            plot_y_z.append(z)
        except IndexError:
            pass 
    print(plot_y_z)
    axs.scatter(plot_x_width, plot_y, 
                # color=cmap(data_cut_idx / len(mass_cuts)),
                # label=r'$10^{%.1f}$'%mass_cuts[cut_idx]+'~'
                # +r'$10^{%.1f}$ '%mass_cuts[cut_idx+1]+'$M_\\odot$/h', 
                s=10)
    
#     tot_plot_x.append(plot_x[np.argsort(plot_x)])
#     tot_plot_y.append(np.array(plot_y)[np.argsort(plot_x)])
# tot_plot_x = [item for sublist in tot_plot_x for item in sublist]
# tot_plot_y = [item for sublist in tot_plot_y for item in sublist]

    print('')

# r_value, p_value = pearsonr(tot_plot_x, tot_plot_y)
# if type == 'DM':
#     reg_color = "blue"
# elif type == 'Hydro':
#     reg_color = "red"
# axs = sns.regplot(x=tot_plot_x, y=tot_plot_y, ci=95, 
#                   scatter=False, line_kws={"color": reg_color})
# plt.title(f'r = {r_value:.2f}, p = {p_value:.3f}', loc='center')

# # Final edit
# axs.set_xlabel(f'Accretion rate width ({acc_width_type})')
# axs.set_ylabel(f"Splashback feature (width)")
# axs.legend()

# # Save the plot
# plt.savefig(f'{root_dir}/{sim}_{type}_sp_width_accret_width_{acc_width_type}')
# plt.close()