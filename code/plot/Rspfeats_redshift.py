from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import numpy as np
import os
import h5py 
import illustris_python as il
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--feat_idx',default=0,type=int) # Feature index [Rsp = 0, depth = 1, width = 2]
args = parser.parse_args()

print('')
print(f'>>> Plot Rsp feats vs redshift <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



boxsize = 205
res = 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)



# Save directory
data_path = f'result/bootstrap_stats/sim_{boxsize}_{res}/Nboots_{args.Nboots}/'
save_dir = f'result/bootstrap_plot/sim_{boxsize}_{res}/Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)



# Inputs
snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]
mass_cut = [11, 11.5, 12, 12.5, 13, 13.5, 14]



feat_idx = args.feat_idx
feats = ['Rsp', 'depth', 'width_dimless', 'width_phys']
print(feats[feat_idx])

# Plot 
fig, axs = plt.subplots(1, 1, dpi=500)

##############################################################################################
### Plot TNG300 ### 
##############################################################################################

# Colour map
cmap = plt.get_cmap('winter', len(mass_cut))

z, all_data = [], []
for snap in snaps: # from high z to low z (present)
    # Load redshifts
    with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z.append(1 / scale_factor - 1)
        
    # Load data
    data = np.load(data_path+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data = data['final_results']
    all_data.append(data)

for cut_idx in range(len(mass_cut)-1): # from low cut to high cut
                     
    # Sort all data
    plot_x, plot_y, plot_y_min, plot_y_max = [], [], [], []
    for snap_idx in range(len(all_data)):
        if all_data[snap_idx].shape[0] > cut_idx:
            plot_y.append(all_data[snap_idx][cut_idx, feat_idx, 1])
            plot_y_min.append(all_data[snap_idx][cut_idx, feat_idx, 0])
            plot_y_max.append(all_data[snap_idx][cut_idx, feat_idx, 2])
            
            plot_x.append(z[snap_idx])
        else:
            pass

    axs.plot(plot_x, plot_y, color=cmap(cut_idx / len(mass_cut)),
             label=r'mass = $10^{%.1f}$'%mass_cut[cut_idx]+'~'+r'$10^{%.1f}$ '%mass_cut[cut_idx+1]+'$M_\\odot$/h')
    axs.fill_between(plot_x, plot_y_min, plot_y_max, color=cmap(cut_idx / len(mass_cut)), alpha=0.2)
    
axs.set_xlabel('z')
if feat_idx == 0:
    axs.set_ylabel(r"$R_{sp}$ [kpc]")
    axs.set_yscale('log')
elif feat_idx == 1:
    axs.set_ylabel("depth")
elif feat_idx == 2:
    axs.set_ylabel(r'width [$R_{200}$]')
elif feat_idx == 3:
    axs.set_ylabel(r"Wdith [kpc]")
    axs.set_yscale('log')

axs.legend()
plt.savefig(os.path.join(save_dir, f'dm_{feats[feat_idx]}_vs_redshift_TNG300'))
plt.close()