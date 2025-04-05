"""We had the mass table but we didn't know 
1. the corresponding halo index at that snapshot;
2. The linking subhalo masses.
This scripts aim to give 
1. the halo index at each snapshot which matches the mass table;
2, the 1st subhalo masses."""

import numpy as np
from tqdm import tqdm
import pandas as pd
import illustris_python as il
import h5py
import argparse
import os

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim_type', default=None, type=str)
args = parser.parse_args()


def get_halo_feat_at_snapX(simpath, df, field):  
    """This function searches the target feature of the subhalos with global indidces.
    e.g. GroupNSub, GroupFristSUb
    """
    work_df = df.copy()
    for irow, row in tqdm(work_df.iterrows(), desc='snaps'):

        snapnum = int(irow[5:])
        print(irow)
        
        # Load subhalo global index 
        GroupData = il.groupcat.loadHalos(simpath+'output', snapnum, fields=field)
        print(f'{field} loaded')
        # -------------------------------------------------------------------------------------
        # Substitute the subhalo global index with the parent halo index in the snapshot X
        work_df.loc[irow] = [GroupData[int(idx)] if pd.notna(idx) else np.nan for idx in row]
        # If any entries are less than 2, replace them with NaN
        work_df.loc[irow] = work_df.loc[irow].where(work_df.loc[irow] >= 2.0, np.nan)
        # -------------------------------------------------------------------------------------
        print(f'{field} assigned')
    return work_df


def get_subhalo_feat_at_snapX(simpath, df, field):  
    """This function searches the target feature of the subhalos with global indidces.
    e.g. SubhaloMass
    """
    for irow, row in tqdm(df.iterrows(), desc='snaps'):

        snapnum = int(irow[5:])
        print(irow)
        
        # Load subhalo global index 
        SubhaloData = il.groupcat.loadSubhalos(simpath+'output', snapnum, fields=field)
        print(f'{field} loaded')
        # -------------------------------------------------------------------------------------
        # Substitute the subhalo global index with the parent halo index in the snapshot X
        df.loc[irow] = [SubhaloData[int(idx)] if pd.notna(idx) else np.nan for idx in row]
        # -------------------------------------------------------------------------------------
        print(f'{field} assigned')
    return df


boxsize, res = 205, 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if args.sim_type == 'DM':
    # DM
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'_DM/'
else:
    # hydro
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'/'


# Load the halp snap index
load_haloIdx_dir = f'result/DMhalo_table_mass/TNG300/sim_{boxsize}_{res}_{args.sim_type}/'
haloIdx_df = pd.read_csv(f'{load_haloIdx_dir}/halo_snap_idx_table.csv', index_col=0)
print("Loaded haloIdx:")
print(haloIdx_df)
print('')


# Convert the snaps to redshift
snaps_list = haloIdx_df.index.tolist()
snaps_dict = {}
for snap in snaps_list:
    with h5py.File(il.snapshot.snapPath(basePath+'/output/', int(snap[5:])), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        snaps_dict[snap] = z


# Load the halo number of subhalos
NSub_df = get_halo_feat_at_snapX(basePath, haloIdx_df, 'GroupNsubs')
print("Loaded GroupNsubs:")
print(NSub_df)
print('')

# Also replace the entries in haloIdx_df where NSub_df is NaN to NaN
haloIdx_df = haloIdx_df.where(pd.notna(NSub_df), np.nan)
print("Loaded haloIdx:")
print(haloIdx_df)
print('')
del NSub_df 


# Load the first Subhalo index
FirstSub_df = get_halo_feat_at_snapX(basePath, haloIdx_df, 'GroupFirstSub')
print("Loaded GroupFirstSub:")
print(FirstSub_df)
print('')
del haloIdx_df

# The index to the second subhalo is simply the first subhalo index + 1
SecondSub_df = FirstSub_df + 1
# Replace the entries where FirstSub_df is NaN with NaN
SecondSub_df = SecondSub_df.where(pd.notna(FirstSub_df), np.nan)
print("Loaded GroupSecondSub:")
print(SecondSub_df)
print('')

# Convert the FirstSub and SecondSub DataFrame to mass
FirstSubMass_df = get_subhalo_feat_at_snapX(basePath, FirstSub_df, 'SubhaloMass')
SecondSubMass_df = get_subhalo_feat_at_snapX(basePath, SecondSub_df, 'SubhaloMass')
print("Loaded SubhaloMass for FirstSub:")
print(FirstSubMass_df)
print("Loaded SubhaloMass for SecondSub:")
print(SecondSubMass_df)
print('')

# Calculate the ratio of subhalo masses
mass_ratio_df = FirstSubMass_df / SecondSubMass_df
del FirstSubMass_df, SecondSubMass_df
print("Calculated mass ratio of FirstSub to SecondSub:")
print(mass_ratio_df)
# If the mass ratio < 7 or > 13, set to nan
mass_ratio_df = mass_ratio_df.where((mass_ratio_df >= 7) & (mass_ratio_df <= 13), np.nan)
print("Filtered mass ratio:")
print(mass_ratio_df)

# Save the data as a 1D array which stores the row index / snapshot, which the
# column entries are most close to 10.
abs_diff = (mass_ratio_df - 10).abs()
select_rows = abs_diff.idxmin(skipna=True)

# Convert the crossing snaps to redshifts
crossing_z = np.array([snaps_dict[snap] if pd.notna(snap) else np.nan
             for snap in select_rows.to_numpy()
             ])
del mass_ratio_df, abs_diff, select_rows

# Check the number of non-NaN crossing redshifts
print(f"Number of non-NaN crossing redshifts: {np.sum(~np.isnan(crossing_z))}")
    
# Load halo data
halos_dir = f'result/DMhalo_density_profiles/TNG300/sim_205_1250_{args.sim_type}/snap_99/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    # Load data
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    data['mergerz'] = crossing_z
    print(data['mergerz'].shape, data['halo_M_Mean200'].shape)
    print(data.keys())
    np.save(halos_dir+fname, data)