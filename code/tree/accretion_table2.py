import illustris_python as il
import h5py
import numpy as np
import argparse
from tqdm import tqdm
from func import calc_tdyn

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim', default='DM-Arepo/MTNG-L500-4320-A/', type=str)
args = parser.parse_args()

print('')
print(f'>>> Accretion rate table 2 <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# BasePath
basePath = f'/virgotng/mpa/MTNG/{args.sim}/'
snaps = [264, 237, 214, 179, 151, 129, 94, 69, 51]

# Load scale factors and redshifts
z_dict, a_dict = {}, {}
for snap in snaps:
    with h5py.File(il.snapshot.snapPath(basePath+'/output/', int(snap)), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        
        z_dict[snap] = z
        a_dict[snap] = scale_factor
    
# Load mass tables
load_dir = f'result/DMhalo_mass_table_new/{args.sim}/'

fname = "snap_{}_FPGrMass.npy"

Nhalos = len(np.load(load_dir+fname.format(snaps[0])))
print(f'Number of halos: {Nhalos}')

accretions = np.zeros((Nhalos, len(snaps) - 1))
for isnap, snap in tqdm(enumerate(snaps[:-1]), desc='calculating accretion rates'):
    
    current_snap_mass = np.load(load_dir+fname.format(snap))
    
    # Find the snapshot which is one t_dyn later
    prev_snap, _ = calc_tdyn(f'{basePath}/output/', z_dict, snap)
    prev_snap_mass = np.load(load_dir+fname.format(prev_snap))
    
    # Get the scalar factors
    a_f, a_i = a_dict[snap], a_dict[prev_snap]
    
    # Compute the accretion rates
    mask = (current_snap_mass != -1) & (prev_snap_mass != -1) & (current_snap_mass != 0) & (prev_snap_mass != 0)
    accretions[mask, isnap] = np.log10(current_snap_mass[mask]/prev_snap_mass[mask]) / np.log10(a_f/a_i)
    print('')
    
np.save(f'{load_dir}/accretion_rates.npy', {'accretions': accretions,
                                            'snaps': snaps[:-1]})