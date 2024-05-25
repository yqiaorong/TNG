import illustris_python as il
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
import argparse
import h5py
from tqdm import tqdm
from func import compt_density_profile_hist, compt_density_profile

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize', default=205, type=int)
parser.add_argument('--res',     default=1250,type=int)
parser.add_argument('--snapnum', default=99,  type=int)
parser.add_argument('--groupnum', default=0,  type=int)
# HaloPos
parser.add_argument('--x', default=0,  type=float) # [ckpc/h]
parser.add_argument('--y', default=0,  type=float) # [ckpc/h]
parser.add_argument('--z', default=0,  type=float) # [ckpc/h]
# Halo mass and radius
parser.add_argument('--M', default=0,  type=float) # [10^10 Msun/h]
parser.add_argument('--R', default=0,  type=float) # [ckpc/h]

parser.add_argument('--method', default=None,  type=str)
parser.add_argument('--save_root_dir', default=None,  type=str)
args = parser.parse_args()

print('')
print(f'>>> DM halo density profile <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Directory where TNG data is stored
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
# Need to pick the box size: 35, 75, or 205 cMpc/h
boxsize = args.boxsize
# Need to pick the resolution level: 540, 1080, or 2160 for the 35 cMpc/h box
#                                    455, 910,  or 1820 for the 75 cMpc/h box
#                                    625, 1250, or 2500 for the 205 cMpc/h box
res = args.res
# Path to the output files for the relevant box size and resolution:
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)

# Specify the snapshot
snapnum = args.snapnum



# Load params
with h5py.File(il.snapshot.snapPath(basePath, snapnum), 'r') as f:
    header = dict(f['Header'].attrs.items())
    scale_factor = header['Time']
    h = header['HubbleParam']
    DMmass = header['MassTable'][1] * 10**10 # [MSun / h]
    # BoxSize = header['BoxSize'] # [ckpc / h]
    


# Select DM halo
groupnum       = args.groupnum
haloPos        = [args.x, args.y, args.z] # [ckpc / h]
halo_M_Mean200 = args.M                   # [10^10 Msun / h]
halo_R_Mean200 = args.R                   # [ckpc / h]
       
# Save data
save_dict = {'halo_R_Mean200': halo_R_Mean200, # [ckpc / h]
             'halo_M_Mean200': halo_M_Mean200, # [10^10 Msun / h]
             'h': h, 'scale_factor': scale_factor, 
             'DMmass': DMmass                  # [MSun / h]
             }




# Set halo spatial boundary
edge = 5 # as a multiple of hal0_R_Mean200
xmin, xmax = haloPos[0] - edge * halo_R_Mean200, haloPos[0] + edge * halo_R_Mean200
ymin, ymax = haloPos[1] - edge * halo_R_Mean200, haloPos[1] + edge * halo_R_Mean200
zmin, zmax = haloPos[2] - edge * halo_R_Mean200, haloPos[2] + edge * halo_R_Mean200



# Load Halos coordinates from snapshot
load_dir = os.path.join(basePath, f'snapdir_{snapnum:03d}')
load_list = os.listdir(load_dir)
load_list = [fname for fname in load_list if fname.endswith('hdf5')]

# Density profiles
densities_bins = []

for idx, file in enumerate(tqdm(load_list)):
    
    # Load coordinates
    snap_path = os.path.join(load_dir, file)
    with h5py.File(snap_path, 'r') as f:

        coords = np.array(f['PartType1/Coordinates']) # [ckpc/h]
        mask = (coords[:,0] >= xmin) & (coords[:,0] <= xmax) & (coords[:,1] >= ymin) & (coords[:,1] <= ymax) & (coords[:,2] >= zmin) & (coords[:,2] <= zmax)
        
        sub_coords = coords[mask]
        del coords, mask
        
    # Compute the density profile
    if sub_coords.shape[0] != 0:
        if args.method == 'hist':
            radial_bins = np.logspace(np.log10(0.01*halo_R_Mean200), 
                                      np.log10(5*halo_R_Mean200), 85) # [ckpc/h]
            DM_masses = [DMmass] * sub_coords.shape[0] # [Msun/h]
            densities = compt_density_profile_hist(sub_coords, DM_masses, haloPos, radial_bins)
        elif args.method == 'old':
            densities, radial_bins = compt_density_profile(sub_coords, haloPos, halo_R_Mean200,
                                                           boxsize*1000)
            densities = densities * DMmass
        densities_bins.append(densities) 
        
densities_bins = np.array(densities_bins)
sum_densities_bins = np.sum(densities_bins, axis=0) # [Msun/h / (ckpc/h)^3]
del densities_bins



# Save directory
save_dir = 'result/'+args.save_root_dir+f'/sim_{boxsize}_{res}/snap_{snapnum}'

if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
save_data_dir = os.path.join(save_dir, 'densities')
if not os.path.exists(save_data_dir):
    os.makedirs(save_data_dir)
    
save_plt_dir = os.path.join(save_dir, 'profiles')
if not os.path.exists(save_plt_dir):
    os.makedirs(save_plt_dir)
    
    
    
# Save data
save_dict['radial_bins'] = radial_bins  # [ckpc/h]
save_dict['densities'] = sum_densities_bins # [(Msun/h)/(ckpc/h)^3]
np.save(os.path.join(save_data_dir, f'halo_{groupnum}'), save_dict)



# # Plot the density profiles    
# plt.figure(figsize=(5,5))
# gs = matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],hspace=0,wspace=0)
# ax = plt.subplot(gs[0])
# ax.plot([halo_R_Mean200, halo_R_Mean200], [0, 10**10], linestyle='--')
# ax.loglog(radial_bins, sum_densities_bins)
# ax.set_xlabel('radius [ckpc/h]')
# ax.set_ylabel(r'density [(M$_{\odot}$/h)/(ckpc/h)$^3$]')
# plt.savefig(os.path.join(save_plt_dir, f'halo_{groupnum}.pdf'))
# plt.close()