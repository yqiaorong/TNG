import illustris_python as il
import h5py
import pandas as pd
import numpy as np
import os
import argparse
from tqdm import tqdm

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim', default='Hydro-Arepo/MTNG-L500-4320-A/', type=str)
args = parser.parse_args()

print('')
print(f'>>> Accretion rate table 2 <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# BasePath
basePath = f'/virgotng/mpa/MTNG/{args.sim}/'
snaps = [264, 214, 151, 94, 69, 51]

# Load scale factors
a = []
for snap in tqdm(snaps, desc='loading scale factors'):
    with h5py.File(il.snapshot.snapPath(basePath+'output', snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
    a.append(scale_factor)
    
# Load mass tables
load_dir = f'result/DMhalo_mass_table_new/{args.sim}/'

fname = "snap_{}_FPGrMass.npy"

Nhalos = len(np.load(load_dir+fname.format(snaps[0])))
print(f'Number of halos: {Nhalos}')

accretions = np.zeros((Nhalos, len(snaps) - 1))
for i in tqdm(range(len(snaps)-1), desc='calculating accretion rates'):
    
    current_snap = np.load(load_dir+fname.format(snaps[i]))
    prev_snap = np.load(load_dir+fname.format(snaps[i+1]))
    
    mask = (current_snap != -1) & (prev_snap != -1) & (current_snap != 0) & (prev_snap != 0)
    accretions[mask, i] = np.log10(current_snap[mask]/prev_snap[mask]) / np.log10(a[i]/a[i+1])

np.save(f'{load_dir}/accretion_rates.npy', accretions)