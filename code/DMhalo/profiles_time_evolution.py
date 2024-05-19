import argparse
import numpy as np
import os
import h5py
from func import float_to_str
import illustris_python as il
from matplotlib import pyplot as plt

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,type=int)
parser.add_argument('--res',    default=1250,type=int)
parser.add_argument('--bin_start', default=3.5, type=float)
parser.add_argument('--bin_end',   default=4, type=float)

parser.add_argument('--root_dir',  default=None, type=str)
args = parser.parse_args()

print('')
print(f'>>> Density profile evolution <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Change bins from floats to strings
bin_start, bin_end = float_to_str(args.bin_start, args.bin_end)

# Load dir
load_dir = f'result/Stacked_{args.root_dir}/sim_{args.boxsize}_{args.res}'

# Walk through the directory
file_list = []
for dirpath, subfolder, fnames in os.walk(load_dir):
    for fname in fnames:
        if fname == f'Bins_{bin_start}_to_{bin_end}.npy':
            file_path = os.path.join(dirpath, fname)
            file_list.append(file_path)
    # Get list of snaps
    if len(subfolder) != 0:
        snaps = subfolder

# Sort index    
snaps = sorted(snaps)
file_list = sorted(file_list)
snap_list = [s for s in snaps if any(s in fpath for fpath in file_list)]

# Set up the final plot
fig, axs = plt.subplots(2, 1, figsize=(10, 15))
axs[0].set_title(f'Mass bin 10^{args.bin_start+10} ~ 10^{args.bin_end+10} (Msun/h) time evolution')

# Set up the colour range
cmap = plt.cm.get_cmap('hsv')
colours = [cmap(i / len(file_list)) for i in range(len(file_list))]



# BasePath
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
boxsize = args.boxsize
res = args.res
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)



for file, snap, c in zip(file_list, snap_list, colours):
    
    # load hubble param and scale factor
    with h5py.File(il.snapshot.snapPath(basePath, snap[-2:]), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        h = header['HubbleParam']

    # Load data
    data = np.load(file, allow_pickle=True).item()
    R200_median = data['R200_median'] # [ckpc/h]
    
    pr = data['profile_radius'] * R200_median * scale_factor / h # [kpc]
    pd = data['profile_densities']

    prf = data['profile_radius_fit'] * R200_median * scale_factor / h # [kpc]
    pdf = data['profile_densities_fit']
    
    sr = data['slopes_radius'] * R200_median * scale_factor / h # [kpc]
    s = data['slopes']

    srf = data['slopes_radius_fit'] * R200_median * scale_factor / h # [kpc]
    sf = data['slopes_fit']
    
    # Plot the density profile
    axs[0].errorbar(pr, pd, fmt='.', color=c, label=f'{snap} data')
    axs[0].errorbar(prf, pdf, color=c, label=f'{snap} fit')
    
    # Plot the fitted gradients
    axs[1].errorbar(sr, s, fmt='.', color=c, label=f'{snap} data')
    axs[1].errorbar(srf, sf, color=c, label=f'{snap} theory')
    
# General settings
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_ylabel(r"$\rho$/$\rho_c$")
axs[0].legend()
# axs[0].set_title(f'Stacked density profiles')

axs[1].set_xscale('log')
axs[1].set_xlabel("r [kpc]")
axs[1].set_ylabel("Slope")
axs[1].legend()
# axs[1].set_title('Finding splashback radius')

plt.tight_layout()
    
# Save directory
save_dir = f'result/Evolution_{args.root_dir}/sim_{args.boxsize}_{args.res}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(os.path.join(save_dir, f'Evolution_{bin_start}_to_{bin_end}'))