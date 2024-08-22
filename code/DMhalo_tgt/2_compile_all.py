import os

snaps = [129, 151, 179, 214, 237, 264]
bin_starts = [[3.5], [3.5, 4], [3.5, 4], [3.5, 4, 4.5], [3.5, 4, 4.5], [3.5, 4, 4.5, 5]]
bin_ends = [[4], [4, 4.5], [4, 4.5], [4, 4.5, 5], [4, 4.5, 5], [4, 4.5, 5, 5.5]]

for isnap, snap in enumerate(snaps):
    num_bins = len(bin_ends[isnap])
    for ibin in range(num_bins):
        os.system(f'python3 code/DMhalo_tgt/compile_chunk_profile.py '+
                  f'--sim DM-Arepo/MTNG-L500-4320-A/output --snapnum {snap} '+
                  f'--bin_start {bin_starts[isnap][ibin]} --bin_end {bin_ends[isnap][ibin]}')
print('finished.')