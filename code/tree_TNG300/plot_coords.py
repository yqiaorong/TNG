import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.cm as cm
from matplotlib import pyplot as plt

print('')
print(f'>>> Plot coordinates of halo lifeline <<<')
print('\nInput arguments:')

# ============================================================================================
# Save directory
# ============================================================================================

# Inpput
DM = 'DM'
load_dir = f'result/DMhalo_coords_table/sim_205_1250_{DM}/'
save_dir = f'{load_dir}/Plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
# ============================================================================================
# Read CM coorinates
# ============================================================================================

df_x = pd.read_csv(load_dir+'coords_0_table.csv', index_col=0)
df_y = pd.read_csv(load_dir+'coords_1_table.csv', index_col=0)
df_z = pd.read_csv(load_dir+'coords_2_table.csv', index_col=0)

# The list of snapshots
snap_list = df_x.index.tolist()
print(snap_list)
num_halos = df_x.shape[1]
print(f'The number of halos: {num_halos}')
column_idx = int(input("Please enter a halo index: "))
halo_global_idx = df_x.columns[column_idx]
print(halo_global_idx)

X = df_x.iloc[:, column_idx]
Y = df_y.iloc[:, column_idx]
Z = df_z.iloc[:, column_idx]

# Get the starting and ending snapshots
for x, snap in zip(X, snap_list):
    if pd.isna(x) == False: 
        start_snap = snap
        break
end_snap = snap_list[-1]

# ============================================================================================
# Plot
# ============================================================================================

# Set the color map
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')

cmap = plt.get_cmap('autumn')
color_list = [cmap(i) for i in np.linspace(1, 0, len(snap_list))]


for x, y, z, snap, color in zip(X, Y, Z, snap_list, color_list):
    if pd.isna(x) == False: 
        ax.scatter(x, y, z, # marker='.', 
                   color=color, label=snap)
ax.set_xlabel("x [ckpc/h]")
ax.set_ylabel("y [ckpc/h]")
ax.set_zlabel("z [ckpc/h]")
ax.set_title(f'Temporal trace of halo {halo_global_idx} from {start_snap} to {end_snap}')
plt.savefig(os.path.join(save_dir, halo_global_idx))  