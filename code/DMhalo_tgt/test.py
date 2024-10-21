import numpy as np

# load_dir = 'result/bootstrap_stats/TNG300/sim_205_1250_Hydro/Nboots_1024/'
# snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]

load_dir = 'result/bootstrap_stats/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
snaps = [129, 151, 179, 214, 237, 264]

for snap in snaps:
    print(snap)
    data = np.load(load_dir + 'snap_%d_Rsp_stats.npy'%snap, allow_pickle=True).item()
    print(data['mass_bins'])
    print('')