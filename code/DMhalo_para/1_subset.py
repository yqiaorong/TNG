import illustris_python as il
import os
import argparse
import numpy as np
from mpi4py import MPI

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',  default=205,  type=int)
parser.add_argument('--res',      default=1250, type=int)
parser.add_argument('--snapnum',  default=99,   type=int)
parser.add_argument('--bin_start',default=1,    type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',  default=3.5,  type=float) # [10^{10+x} Msun/h]
parser.add_argument('--method', default='old',  type=str)
parser.add_argument('--save_root_dir',default='DMhalo_density_profiles',type=str)
args = parser.parse_args()

print('')
print(f'>>> DM halo density profiles subset <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Specify the snapshot
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(args.boxsize,args.res)
snapnum = args.snapnum
boxsize = args.boxsize
res     = args.res



# Save root dir
save_root_dir = args.save_root_dir +'_'+args.method
    
    

# Select a subset of DM halos from groupcat
group_fields = ['GroupPos', 'Group_M_Mean200', 'Group_R_Mean200']
Halos = il.groupcat.loadHalos(basePath, snapnum, fields=group_fields)

GroupPos        = Halos['GroupPos']        # [ckpc/h]
Group_M_Mean200 = Halos['Group_M_Mean200'] # [10^10 MSun/h]
Group_R_Mean200 = Halos['Group_R_Mean200'] # [ckpc/h]

# Using physical mass to select subset
subset_idx = np.where((Group_M_Mean200 >= 10**args.bin_start) & (Group_M_Mean200 < 10**args.bin_end))[0]
Ngroups_subset = subset_idx.shape[0]
print(f'In total, {Ngroups_subset} DM halos with mass 10^{args.bin_start+10} ~ 10^{args.bin_end+10} MSun in at snap {snapnum}')



# Parallel computing
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(f'size: {size}, rank: {rank}')

# Define initial data (the data you want to split) on the root MPI process
if rank == 0:
    data = np.linspace(0, Ngroups_subset, size+1).astype(int)[:-1]
else:
    data = None
print(f'scatter data: {data}')

# Scatter data to all MPI processes
start_idx_per_core = comm.scatter(data, root=0)

# Debugging: print the indices received by each process
print(f'Process {rank} received start_idx_per_core: {start_idx_per_core}')

# Apply the computation
comm.Barrier()



# Iterate over DM halos
for idx in subset_idx[start_idx_per_core:]:
    if not os.path.exists(f'result/{save_root_dir}/sim_{boxsize}_{res}/snap_{snapnum}/densities/halo_{idx}.npy'):
        # Round values 
        x, y, z = np.round(GroupPos[idx, 0].item(), 0), np.round(GroupPos[idx, 1].item(), 0), np.round(GroupPos[idx, 2].item(), 0)
        R = np.round(Group_R_Mean200[idx].item(), 0)
        # Run the script
        os.system(f'python3 code/DMhalo/one_halo_hist.py'+
            f' --boxsize {boxsize} --res {res} --snapnum {snapnum} --groupnum {idx}'+
            f' --x {x} --y {y} --z {z}'+
            f' --M {Group_M_Mean200[idx]} --R {R}'+
            f' --save_root_dir {save_root_dir} --method {args.method}')
    else:
        print(f'Processor {rank}: At snap {snapnum}, DM halo local index {idx}/{Ngroups_subset-1} already exists.')

print(f'All DM halos in the subset at snap {snapnum} are finished.')

# Gather results on the root process