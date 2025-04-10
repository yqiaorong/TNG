"""This script creats an identical halo table of halo's pH and find halo formation time using pH."""

import numpy as np
import pandas as pd
import illustris_python as il
import h5py
import os
import argparse
from colossus.cosmology import cosmology
from colossus.lss import peaks

cosmology.setCosmology('planck15')


# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim', default='TNG300/sim_205_1250_Hydro/', type=str) # [TNG300/MTNG]
args = parser.parse_args()

print('')
print(f'>>> Assign peak height to halo history <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if 'DM' in args.sim:
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'_DM/output/'
elif 'Hydro' in args.sim:
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'/output/'
    

# Load the halo mass table
load_halo_dir = f'result/DMhalo_table_mass/{args.sim}/'
haloMass_df = pd.read_csv(f'{load_halo_dir}/halo_mass_table.csv', index_col=0) # Comoving mass
print("Loaded haloMass:")
print(haloMass_df)
print('')        


# Create an empty dataframe as haloMass_df
pH_df = pd.DataFrame(index=haloMass_df.index, columns=haloMass_df.columns)


# Load halo data
snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]
snaps_dict = {}
for snap in snaps:
    
    # Load redshift
    with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        h = header['HubbleParam']
        snaps_dict[f'snap_{snap}'] = z
        
    # Convert the halo mass from comoving to physical!!!
    Massp = haloMass_df.loc[f'snap_{snap}'].values / h # [10^10 MSun]
    # Caculate the peak height
    peakHeight = peaks.peakHeight(Massp*10**10, z)
    # Assign pH to df
    pH_df.loc[f'snap_{snap}'] = peakHeight
print("halo peakHeight:")            
print(pH_df)
print('')


# Save the data as a 1D array which stores the row index / snapshot, which the
# column entries are most close to 1.
abs_diff = (pH_df - 1).abs()
select_rows = abs_diff.idxmin(skipna=True)

# Convert the crossing snaps to redshifts
crossing_z = np.array([snaps_dict[snap] if pd.notna(snap) else np.nan
             for snap in select_rows.to_numpy()
             ])
del abs_diff, select_rows

# Check the number of non-NaN crossing redshifts
print(f"Number of non-NaN crossing redshifts: {np.sum(~np.isnan(crossing_z))}")
    
# Load halo data
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_99/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    # Load data
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    data['formz'] = crossing_z
    print(data['formz'].shape, data['halo_M_Mean200'].shape)
    data.pop('formation_time', None)
    print(data.keys())
    np.save(halos_dir+fname, data)
    print('data saved! ')