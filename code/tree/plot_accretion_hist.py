import pandas as pd
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import h5py
import illustris_python as il
import numpy as np



DM = ''
boxsize, res = 205, 1250

data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)

load_dir = f'result/DMhalo_mass_table/sim_{boxsize}_{res}{DM}/'

mass_df = pd.read_csv(load_dir+'mass_table.csv', index_col=0)
accretion_df = pd.read_csv(load_dir+'accretion_table.csv', index_col=0)

snap_list = mass_df.index

# Load redshift
redshifts = []
for snap in snap_list:
    with h5py.File(il.snapshot.snapPath(basePath, int(snap[5:])), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        redshifts.append(z)

for idx in range(len(snap_list)):
    if idx != 0:
        rate = accretion_df.iloc[idx, :]
        rate = rate.replace([np.inf, -np.inf], np.nan)
        rate = rate.dropna().to_numpy()
        median = np.median(rate)

        plt.figure()
        hist = plt.hist(rate, label=f'z = {np.round(redshifts[idx], 2)}')[0]
        plt.plot([median, median], [0, max(hist)], linestyle='dashed', color='red', 
                  label=f'median = {np.round(median, 3)}')
        plt.legend(loc='best')
        plt.savefig(load_dir+f'snap_{snap_list[idx][5:]}')
    