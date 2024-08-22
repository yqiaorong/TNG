from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import numpy as np
import os
import h5py 
import illustris_python as il
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default=None,type=str)
parser.add_argument('--feat_idx',default=0,   type=int) # Feature index [Rsp = 0, depth = 1, width = 2]
args = parser.parse_args()

print('')
print(f'>>> Plot Rsp feats vs redshift <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

basePath = f'/virgotng/mpa/MTNG/{args.sim}'

# Save directory
load_dir = f'result/bootstrap_stats/{args.sim}'
save_dir = f'result/bootstrap_plot/{args.sim}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Inputs
snaps = [129, 151, 179, 214, 237, 264]
mass_cut = [13.5, 14, 14.5, 15]

feat_idx = args.feat_idx
feats = ['Rsp', 'depth', 'width-dimless', 'width-phys']
factors = [1000, 1, 1, 1000]
print(feats[feat_idx])

# Plot 
fig, axs = plt.subplots(1, 1, dpi=500)
# Colour map
cmap = plt.get_cmap('winter', len(snaps))

z, all_data = [], []
for snap in snaps:
    # Load redshifts
    with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z.append(1 / scale_factor - 1)
        
    # Load data
    data = np.load(load_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data = data['final_results']
    all_data.append(data)

for cut_idx in range(len(mass_cut)-1): # from low cut to high cut
    
    # Because of difference in unit, Rsp should be multipled by 1000 from [Mpc] to [kpc]
    factor = factors[feat_idx]
        
    # Sort all data
    plot_x, plot_y, plot_y_min, plot_y_max = [], [], [], []
    for snap_idx in range(len(all_data)):
        if all_data[snap_idx].shape[0] > cut_idx:
            plot_y.append(factor*all_data[snap_idx][cut_idx, feat_idx, 1])
            plot_y_min.append(factor*all_data[snap_idx][cut_idx, feat_idx, 0])
            plot_y_max.append(factor*all_data[snap_idx][cut_idx, feat_idx, 2])
            
            plot_x.append(z[snap_idx])
        else:
            pass
    
    # axs.errorbar(plot_x, plot_y, 
    #              yerr = [np.abs(plot_y_min-plot_y), np.abs(plot_y_max-plot_y)],
    #              label=r'mass = $10^{%.1f}$'%mass_cut[cut_idx]+'~'+r'$10^{%.1f}$ '%mass_cut[cut_idx+1]+f'$M_\\odot$/h)
    axs.plot(plot_x, plot_y, color=cmap(cut_idx / len(mass_cut)),
             label=r'mass = $10^{%.1f}$'%mass_cut[cut_idx]+'~'+r'$10^{%.1f}$ '%mass_cut[cut_idx+1]+'$M_\\odot$/h')
    axs.fill_between(plot_x, plot_y_min, plot_y_max, 
                     color=cmap(cut_idx / len(mass_cut)), alpha=0.2)
    
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

if args.sim.startswith('DM-Arepo'):
    plot_name = f'dm_{feats[feat_idx]}_vs_redshift_MTNG-DM'
elif args.sim.startswith('Hydro-Arepo'):
    plot_name = f'dm_{feats[feat_idx]}_vs_redshift_MTNG-Hydro'    
plt.savefig(f'{save_dir}/{plot_name}')
plt.close()
