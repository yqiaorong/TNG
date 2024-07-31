from matplotlib import pyplot as plt 
plt.style.use('code/plot/style.mplstyle')
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
data_path = f'result/bootstrap'
save_dir = f'{data_path}/sim_{boxsize}_{res}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


# Inputs
snaps = [33, 40, 50, 67, 78, 99]
mass_cut = [11, 11.5, 12, 12.5, 13, 13.5, 14]



feats_idx = args.feat_idx
feats = ['Rsp', 'depth', 'width']
print(feats[feats_idx])

# Plot 
fig, axs = plt.subplots(1, 1, dpi=500)

z, all_data = [], []
for i, snap in enumerate(snaps):
    # Load redshifts
    with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z.append(1 / scale_factor - 1)
        
    # Load data
    data = np.load(data_path+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data = data ['final_results']
    all_data.append(data)
all_data = np.array(all_data)

for i in range(len(mass_cut)-1):
    # axs.errorbar(z, all_data[:, 1, i, 0], 
    #              yerr = [np.abs(all_data[:, 0, i, 0]-all_data[:, 1, i, 0]), 
    #                      np.abs(all_data[:, 2, i, 0]-all_data[:, 1, i, 0])],
    #              color=colours[i], label=f'mass = 10^{mass_cut[i]} MSun/h')
    axs.plot(z, all_data[:, i, feats_idx, 1], 
             label=r'mass = $10^{%.1f}$'%mass_cut[i]+'~'+r'$10^{%.1f}$ '%mass_cut[i+1]+f'$M_\\odot$/h')
    axs.fill_between(z, all_data[:, i, feats_idx, 0], all_data[:, i, feats_idx, 2], alpha=0.2)
    
axs.set_xlabel('z')
if args.feat_idx == 0:
    axs.set_ylabel(r"$R_{sp}$ [kpc]")
else:
    axs.set_ylabel(feats[feats_idx])
if args.feat_idx == 0:
    axs.set_yscale('log')
axs.legend()
# plt.tight_layout() # incompatible with pltstyle
plt.savefig(os.path.join(save_dir, f'dm_{feats[feats_idx]}_vs_redshift_TNG300'))
plt.close()