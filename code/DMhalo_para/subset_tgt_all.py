import os
import argparse
from mpi4py import MPI
import numpy as np


# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',          default=None,                     type=str)
parser.add_argument('--snapnum',      default=264,                      type=int)
parser.add_argument('--bin_start',    default=None,                     type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',      default=None,                     type=float) # [10^{10+x} Msun/h]
parser.add_argument('--save_root_dir',default='DMhalo_density_profiles',type=str)
args = parser.parse_args()

# Set up directory
basePath = f'/nfs/mvogelsblab001/Users/s_qyu/MTNG/{args.sim}/'
snapnum = args.snapnum

load_dir = os.path.join(basePath, f'snapdir_{snapnum:03d}')
load_list = os.listdir(load_dir)
load_list = [fname for fname in load_list if fname.startswith(f'snapshot_{snapnum:03d}') and fname.endswith('hdf5')]
print(len(load_list))



# Parallel computing
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(f'size: {size}, rank: {rank}')



# Define initial data (the data you want to split) on the root MPI process
if rank == 0:
    # data = list(range(0, len(load_list), int(math.floor(len(load_list)/size))))
    data = np.linspace(0, len(load_list), size+1)
    data = np.round(data).astype(int).tolist()
    data = data[:size]
else:
    data = None
print(f'scatter data: {data}')

# Scatter data to all MPI processes
start_idx_per_core = comm.scatter(data, root=0)

# Debugging: print the indices received by each process
print(f'Process {rank} received start_idx_per_core: {start_idx_per_core}')

# Apply the computation
comm.Barrier()



# Running
for chunk in range(start_idx_per_core, start_idx_per_core+int(len(load_list)/size)):
    if not os.path.exists(f'result/{args.save_root_dir}/{args.sim}/snap_{snapnum}/intermediate_densities/'+
                          f'bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}/chunk-{chunk}.npy'):
        os.system('python3.11 code/DMhalo_tgt/subset_tgt_chunk-gpu.py'+
             f' --sim {args.sim} --snapnum {args.snapnum} --chunk_idx {chunk}'+
             f' --bin_start {args.bin_start} --bin_end {args.bin_end}'
             f' --save_root_dir {args.save_root_dir}')
    else:
        print(f'chunk-{chunk}_bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}.npy: exist')
print(f'All files in the chunk {chunk} in snap {args.snapnum} are finished.')