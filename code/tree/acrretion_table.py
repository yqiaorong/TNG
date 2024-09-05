import illustris_python as il
import h5py
import pandas as pd
import numpy as np
import os

# input
DM = ''
boxsize = 205
res = 1250

# BasePath
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+f'{DM}/output/'



# Load mass csv
df = pd.read_csv(f'result/DMhalo_mass_table/sim_{boxsize}_{res}{DM}/mass_table.csv', index_col=0)

# Load scale factors
snap_list = df.index
print(snap_list)
a = []
for snap in snap_list:
    # Load scale factors
    with h5py.File(il.snapshot.snapPath(basePath, int(snap[5:])), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
    a.append(scale_factor)

# Create accretion rate df
rate_df = pd.DataFrame(index=df.index, columns=df.columns)

# Compute accretion rates
for isnap in range(len(snap_list)-1):
    
    mass_f, mass_i = df.iloc[isnap+1], df.iloc[isnap]
    a_f, a_i = a[isnap+1], a[isnap]
    
    mask = mass_f.notna() & mass_i.notna()
    
    # Calculate the accretion rate
    rate = np.log10(mass_f[mask]/mass_i[mask]) / np.log10(a_f/a_i)
    
    # Assign the calculated rate back to the new DataFrame
    
    rate_df.loc[snap_list[isnap+1], mask.index[mask]] = rate

# Save the accretion rate df
save_mass_dir = f'result/DMhalo_mass_table/sim_{boxsize}_{res}{DM}'
if not os.path.exists(save_mass_dir):
    os.makedirs(save_mass_dir)
rate_df.to_csv(f'{save_mass_dir}/accretion_table.csv', index=snap_list)