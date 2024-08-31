import os
import numpy as np

snaps = [8, 17, 21, 25, 33, 40, 50, 67, 78, 99]
mass_bins = np.arange(1, 4.5, 0.5)

for snap in snaps:
    for bin_start in mass_bins:
        if not os.path.exists(f'result/DMhalo_density_profiles/TNG300/sim_205_1250/snap_{snap}/final_densities/bin-{int(bin_start*10)}-{int(bin_start*10+5)}.npy'):
            os.system(f'python3 code/DMhalo/compile_chunk_profile.py --snapnum {snap} --bin_start {bin_start} --bin_end {bin_start+0.5}')