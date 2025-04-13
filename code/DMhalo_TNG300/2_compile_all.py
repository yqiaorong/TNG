"""All the data computed in this script is saved in result/DMhalo_density_profiles_raw2"""

import os
import numpy as np

snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]

# Input 
DM = 'Hydro' # [Hydro / DM]

for snap in zip(snaps):
    if not os.path.exists(f'result/DMhalo_density_profiles_raw2/TNG300/sim_205_1250_{DM}/snap_{snap}/final_densities/bin-10-50.npy'):
        os.system(f'python3 code/DMhalo/compile_chunk_profile.py --snapnum {snap} --DM {DM}')