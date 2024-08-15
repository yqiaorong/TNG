"""Add data to existing data files."""

import os
import h5py
import numpy as np
import illustris_python as il
from tqdm import tqdm

odir = f'result/bootstrap_stats/sim_205_1250/'

boxsize = 205
res = 1250
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/' + 'L%dn%dTNG/output'%(boxsize,res)

snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]
for snap in tqdm(snaps):
    
    # Load redshift values
    with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        Parameters = dict(f['Parameters'].attrs.items())

        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        h = Parameters['HubbleParam']
        
    data = np.load(f'{odir}snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data['z'] = z
    data['h'] = h
    np.save(f'{odir}/snap_{snap}_Rsp_stats.npy', data)