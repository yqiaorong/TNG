import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import os

# Specify the snapshot
basePath = '/nfs/mvogelsblab002/Users/s_qyu/TNG300-1/output'
snapNum = 99

# Load Halos from groupcat
Group_M_Mean200 = il.groupcat.loadHalos(basePath, snapNum, fields='Group_M_Mean200')

# Save directory
save_dir = 'result/DM halos mass histogram'
if os.path.isdir(save_dir) == False:
    os.makedirs(save_dir)
    
# Histogram of halos mass
plt.figure()
plt.hist(Group_M_Mean200, bins=np.logspace(0.01, 6, 50))
plt.xlabel('Mass [10^10 MSun / h]')
plt.ylabel('Frequency')
plt.yscale('log')
plt.xscale('log')
plt.title(f"DM halos' mass histogram at snap {snapNum}")
plt.savefig(f'result/DM halos mass histogram/snap_{snapNum}')