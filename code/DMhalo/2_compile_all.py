"""All the data computed in this script is saved in result/DMhalo_density_profiles_raw2"""

import os
import numpy as np

snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]
bin_ends = [2, 3, 3, 3.5, 3.5, 4, 4, 4, 4.5, 4.5, 4.5]

# Input 
DM = 'DM' # [Hydro / DM]

for snap, bin_end in zip(snaps, bin_ends):
    if not os.path.exists(f'result/DMhalo_density_profiles/TNG300/sim_205_1250_{DM}/snap_{snap}/final_densities/bin-10-{int(bin_end*10+5)}.npy'):
        os.system(f'python3 code/DMhalo/compile_chunk_profile.py --snapnum {snap} --DM {DM} --bin_start 1 --bin_end {bin_end}')