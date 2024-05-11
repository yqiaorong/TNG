import illustris_python as il
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import argparse
import h5py
from func import compt_density_profile

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize', default=205, type=int)
parser.add_argument('--res',     default=1250,type=int)
parser.add_argument('--snapnum', default=99,  type=int)
parser.add_argument('--groupnum',default=0,   type=int)
args = parser.parse_args()

print('')
print(f'>>> DM halo density profile <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Directory where TNG data is stored
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
# Need to pick the box size: 35, 75, or 205 Mpc/h
boxsize = args.boxsize
# Need to pick the resolution level: 540, 1080, or 2160 for the 35Mpc/h box
#                                    455, 910,  or 1820 for the 75Mpc/h box
#                                    625, 1250, or 2500 for the 205Mpc/h box
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
    DMmass = header['MassTable'][1] * 10**10 # unit MSun



# Load Halos from groupcat
group_fields = ['GroupPos', 'Group_M_Mean200', 'Group_R_Mean200']
Halos = il.groupcat.loadHalos(basePath, snapnum, fields=group_fields)

GroupPos        = Halos['GroupPos']        / h * scale_factor
Group_M_Mean200 = Halos['Group_M_Mean200'] / h 
Group_R_Mean200 = Halos['Group_R_Mean200'] / h * scale_factor
del Halos

# Select DM halo
groupnum       = args.groupnum
haloPos        = GroupPos[groupnum]
halo_M_Mean200 = Group_M_Mean200[groupnum]
halo_R_Mean200 = Group_R_Mean200[groupnum]
del GroupPos, Group_M_Mean200, Group_R_Mean200



# Save data
save_dict = {'halo_R_Mean200': halo_R_Mean200, 'halo_M_Mean200': halo_M_Mean200,
             'h': h, 'scale_factor': scale_factor, 'DMmass': DMmass}

# Set halo spatial boundary
edge = 5 # as a multiple of hal0_R_Mean200
xmin, xmax = haloPos[0] - edge * halo_R_Mean200, haloPos[0] + edge * halo_R_Mean200
ymin, ymax = haloPos[1] - edge * halo_R_Mean200, haloPos[1] + edge * halo_R_Mean200
zmin, zmax = haloPos[2] - edge * halo_R_Mean200, haloPos[2] + edge * halo_R_Mean200

# Load Halos coordinates from snapshot
# Coordinates = il.snapshot.loadSubset(basePath, snapNum, 'dm', ['Coordinates'], float32=True)

load_dir = os.path.join(basePath, f'snapdir_{snapnum:03d}')
load_list = os.listdir(load_dir)
load_list = [fname for fname in load_list if fname.endswith('hdf5')]

rho_bins, r_bins = [], []
for idx, file in enumerate(load_list):
    
    # Load coordinates
    snap_path = os.path.join(load_dir, file)
    with h5py.File(snap_path, 'r') as f:

        coords = np.array(f['PartType1/Coordinates']) / h * scale_factor
        mask = (coords[:,0] >= xmin) & (coords[:,0] <= xmax) & (coords[:,1] >= ymin) & (coords[:,1] <= ymax) & (coords[:,2] >= zmin) & (coords[:,2] <= zmax)
        
        sub_coords = coords[mask]
        del coords, mask
    
    # Compute the density profile
    if sub_coords.shape[0] != 0:
        rho, r = compt_density_profile(sub_coords, haloPos, halo_R_Mean200)
        rho = rho * DMmass
        rho_bins.append(rho)
        r_bins = r
    del sub_coords
rho_bins = np.array(rho_bins)
# Get the final density profile
densities = np.sum(rho_bins, axis=0)



# Save directory
save_dir = f'result/DMhalo_density_profiles/sim_{args.boxsize}_{args.res}/snap_{snapnum}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
save_data_dir = os.path.join(save_dir, 'densities')
if not os.path.exists(save_data_dir):
    os.makedirs(save_data_dir)
    
save_plt_dir = os.path.join(save_dir, 'profiles')
if not os.path.exists(save_plt_dir):
    os.makedirs(save_plt_dir)
    
# Save data
save_dict['radial_bins'] = r_bins
save_dict['densities'] = densities
np.save(os.path.join(save_data_dir, f'halo_{groupnum}'), save_dict)

# Plot the density profiles    
plt.figure(figsize=(5,5))
gs = matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],hspace=0,wspace=0)
ax = plt.subplot(gs[0])
ax.plot([halo_R_Mean200, halo_R_Mean200], [0, 10**10], linestyle='--')
ax.loglog(r_bins, densities)
ax.set_xlabel('radius [kpc]')
ax.set_ylabel(r'density [M$_{\odot}$/kpc$^3$]')
plt.savefig(os.path.join(save_plt_dir, f'halo_{groupnum}.pdf'))
plt.close()