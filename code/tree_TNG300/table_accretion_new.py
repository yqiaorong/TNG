import illustris_python as il
import h5py
import pandas as pd
import numpy as np
import os
import argparse
from func import calc_tdyn

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim_type', default=None, type=str)
args = parser.parse_args()

# input
sim_type = args.sim_type
boxsize = 205
res = 1250

# BasePath
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if sim_type == 'DM':
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'_DM/output/'
elif sim_type == 'Hydro':
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'/output/'



# Load mass csv
df = pd.read_csv(f'result/DMhalo_table_mass/TNG300/sim_{boxsize}_{res}_{sim_type}/halo_mass_table.csv', index_col=0)

# Create a dictionary of redshifts
snap_list = df.index
z_dict, a_dict = {}, {}
for snap in snap_list:
    with h5py.File(il.snapshot.snapPath(basePath, int(snap[5:])), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        
        z_dict[snap] = z
        a_dict[snap] = scale_factor
        
        
        
# Create accretion rate df
rate_df = pd.DataFrame(index=df.index, columns=df.columns)

for snap in snap_list:
    print('current snapshot: ', snap)
    
    # Find the snapshot which is one t_dyn later
    prev_snap, _ = calc_tdyn(basePath, z_dict, int(snap[5:]))
    
    if prev_snap != snap:
        # Load the corresponding mass
        mass_f = df.loc[snap]
        mass_i = df.loc[prev_snap]
    
        # Calculate the accretion rate
        a_f, a_i = a_dict[snap], a_dict[prev_snap]
        print(a_f, a_i)
        # Mask out the nan values and zero values
        mask = (mass_f.notna() & mass_i.notna()) & (mass_i != 0) & (mass_f != 0)
        rate = np.log10(mass_f[mask]/mass_i[mask]) / np.log10(a_f/a_i)
        # Check if rate contains inf
        if np.isinf(rate).any():
            print('Inf values detected')
        
        # Assign the calculated rate back to the new DataFrame
        rate_df.loc[snap, mask.index[mask]] = rate
    print('')
    
rate_df.index = snap_list
    
    

# Save the accretion rate df
save_mass_dir = f'result/DMhalo_table_mass/TNG300/sim_{boxsize}_{res}_{DM}'
if not os.path.exists(save_mass_dir):
    os.makedirs(save_mass_dir)
rate_df.to_csv(f'{save_mass_dir}/accretion_table_new.csv', index=True)