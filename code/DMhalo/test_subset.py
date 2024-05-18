import matplotlib.pyplot as plt
import illustris_python as il
import argparse
import numpy as np
import os
import h5py
from unyt import Msun, g, cm, kiloparsec

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,   type=int)
parser.add_argument('--res',     default=2500, type=int)
parser.add_argument('--snapnum', default=99,   type=int)
parser.add_argument('--mass_range',default=3.5,type=float) # [10^{10+x} Msun]
args = parser.parse_args()

print('')
print(f'>>> Test DM halo density profiles subset <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

snapnum = args.snapnum

# Directory where TNG data is stored
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
# Need to pick the box size: 35, 75, or 205 Mpc/h
boxsize = args.boxsize
# Need to pick the resolution level: 540, 1080, or 2160 for the 35Mpc/h box
#                                    455, 910,  or 1820 for the 75Mpc/h box
#                                    625, 1250, or 2500 for the 205Mpc/h box
res = args.res
# Path to the output files for the relevant box size and resolution:
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)



# Load params
with h5py.File(il.snapshot.snapPath(basePath, snapnum), 'r') as f:
    header = dict(f['Header'].attrs.items())
    scale_factor = header['Time']
    h = header['HubbleParam'] # unit [100 * km / megaparsec / second]
    BoxSize = header['BoxSize'] # [ckpc/h]

    

with h5py.File(f'data/halo_data_res{args.res}_snap{args.snapnum}.hdf5', 'r') as f:
    HaloIndices = f['HaloIndices']
    
    # Raw mass table
    Group_M_Mean200 = f['HaloM200mean'] * g / (10**10 * Msun) * h # unit [10^10 MSun/h]
    subset_idx = [idx for idx, mass in enumerate(Group_M_Mean200)
                if mass >= 10**args.mass_range and mass < 10**(args.mass_range+0.5)]
    print(f'Number of halos in the subset: {len(subset_idx)}')

    # Scaled table
    Group_R_Mean200 = f['HaloR200mean'] * cm / kiloparsec * h / scale_factor # unit [ckpc/h]
    GroupPos        = f['HaloPositions'] * cm / kiloparsec * h / scale_factor  # unit [ckpc/h]

    
    # Iterate over DM halos
    for idx in subset_idx:
        path = f'result/test_res{res}_snap{snapnum}/'+f'sim_{args.boxsize}_{res}/snap_{snapnum}/densities/halo_{HaloIndices[idx]}.npy'
        if not os.path.exists(path):
            # Round values 
            x, y, z = np.round(GroupPos[idx, 0].item(), 0), np.round(GroupPos[idx, 1].item(), 0), np.round(GroupPos[idx, 2].item(), 0)
            R = np.round(Group_R_Mean200[idx].item(), 0)
            os.system(f'python3 code/DMhalo/one_halo_test.py'+
                f' --boxsize {args.boxsize} --res {res} --snapnum {snapnum}'+
                f' --groupnum {HaloIndices[idx]}'+
                f' --x {x} --y {y} --z {z}'+
                f' --M {Group_M_Mean200[idx].item()} --R {R}'+
                f' --save_root_dir /test_res{res}_snap{snapnum}/')

    print(f'All DM halos in the subset at snap {snapnum} are finished.')
    
    
    
# Histogram of raw halos mass
plt.figure()
hist_values = plt.hist(Group_M_Mean200, bins=np.logspace(3, 5, 50))[0]
plt.xlabel('Mass [10^10 MSun]')
plt.ylabel('Frequency')
plt.yscale('log')
plt.xscale('log')
plt.title("DM halos' mass histogram")
plt.savefig('result/DM halos mass histogram/test_mass_hist')