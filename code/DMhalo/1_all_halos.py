import illustris_python as il
import os

# Specify the snapshot
basePath = '/nfs/mvogelsblab002/Users/s_qyu/TNG300-1/output'
snapNum = 99

# Get the total number of DM halos
Header = il.groupcat.loadHeader(basePath, snapNum)
Ngroups_Total = Header['Ngroups_Total']

# Iterate over DM halos
for halo_idx in range(Ngroups_Total):
    if not os.path.exists(f'result/DM_halo_density_profiles/snap_{snapNum}/halo_{halo_idx}.npy'):
        os.system(f'python3 code/DMhalo/one_halo.py --halo_idx {halo_idx}')
    else:
        print(f'At snap {snapNum}, DM halo local index {halo_idx+1}/{Ngroups_Total} already exists.')

print(f'All halos at snap {snapNum} are finished.')