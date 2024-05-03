# %%
### check chunk file number in one snapshot ###

import os

output_dir = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L205n1250TNG/output'

# choose snapshot
snap = 78

snap_dir = os.path.join(output_dir, f'snapdir_{snap:03d}')
# Check the number of hdf5 chunk file
snap_list = os.listdir(snap_dir)
snap_list = sorted(snap_list)

num_sfile = 0
for file in snap_list:
    if file.endswith('hdf5'):
        num_sfile += 1
print(f'The number of snap chunk file in snap {snap}: {num_sfile}/600')

group_dir = os.path.join(output_dir, f'groups_{snap:03d}')
# Check the number of hdf5 chunk file
group_list = os.listdir(group_dir)
group_list = sorted(group_list)

num_gfile = 0 
for file in group_list:
    if file.endswith('hdf5'):
        num_gfile += 1
print(f'The number of group chunk file in snap {snap}: {num_gfile}/600')

# %%
### Rename files

fdir = os.path.join(output_dir, 'snapdir_099')
flist = os.listdir(fdir)
for f in flist:
    if f.startswith('snapshot-99'):
        new_f = f.replace('snapshot-99', 'snap_099')
        old_path = os.path.join(fdir, f)
        new_path = os.path.join(fdir, new_f)
        os.rename(old_path, new_path)
del fdir, flist
        
# %%
### Check number of DM halo computed
import os

parent_dir = '/n/home01/sqyu'
halo_dir = os.path.join(parent_dir, 'AstroLab/TNG300-1/result/DM_halo_density_profiles',
                        'snap_99')
halo_list = os.listdir(halo_dir)
print(len(halo_list))
del halo_dir, halo_list
# %%
