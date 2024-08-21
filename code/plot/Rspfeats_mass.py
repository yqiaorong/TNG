from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import numpy as np
import os
import h5py 
import illustris_python as il
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--Nboots',  default=None,type=int)
parser.add_argument('--feat_idx',default=0,   type=int) # Feature index [Rsp = 0, depth = 1, width = 2]
args = parser.parse_args()

print('')
print(f'>>> Plot Rsp feats vs mass <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

boxsize = 205
res = 1250
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/' + 'L%dn%dTNG/output'%(boxsize,res)



# Save directory
data_path = f'result/bootstrap_stats/sim_{boxsize}_{res}/Nboots_{args.Nboots}/'
save_dir = f'result/bootstrap_plot/sim_{boxsize}_{res}/Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


# Input
snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]
mass_cut = [10**11, 10**11.5, 10**12, 10**12.5, 10**13, 10**13.5]


feat_idx = args.feat_idx
feats = ['Rsp', 'depth', 'width_dimless', 'width_phys']
print(feats[feat_idx])

# Plot 
fig, axs = plt.subplots(1, 1, 
                        # figsize=(10, 8), 
                        dpi=500)
# Colour map
cmap = plt.get_cmap('winter', len(snaps))

z = []
for isnap, snap in enumerate(snaps): # from high z to low z
    # Load redshifts
    with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z.append(1 / scale_factor - 1)
        
    # Load data
    data = np.load(data_path+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data = data['final_results']
    num_cut = data.shape[0]
    
    factors = [1000, 1, 1, 1000]
    data[:, feat_idx, :] = data[:, feat_idx, :]*factors[feat_idx] # convert Rsp in [Mpc] to [kpc]
    
    # axs.errorbar(mass_cut[:num_cut], data[:, feats_idx, 1], 
    #              yerr = [np.abs(data[:, feats_idx, 1]-data[:, feats_idx, 0]), 
    #                      np.abs(data[:, feats_idx, 1]-data[:, feats_idx, 2])],
    #              label=f'z = {np.round(z[i], 3)}')
    axs.plot(mass_cut[:num_cut], data[:, feat_idx, 1], color=cmap(isnap / len(snaps)),
             label=f'z = {np.round(z[isnap], 1)}')
    axs.fill_between(mass_cut[:num_cut], data[:, feat_idx, 0], data[:, feat_idx, 2], 
                     color=cmap(isnap / len(snaps)), alpha=0.2)
    
axs.set_xscale('log')
axs.set_xlabel('Mass [$M_\\odot$/h]')
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
# plt.tight_layout() # incompatible with pltstyle
plt.savefig(f'{save_dir}/dm_{feats[feat_idx]}_vs_mass_TNG300')
plt.close()