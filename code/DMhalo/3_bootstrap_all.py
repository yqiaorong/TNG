import os
from mpi4py import MPI

# Parallel computing
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(f'size: {size}, rank: {rank}')

# Define initial data (the data you want to split) on the root MPI process
if rank == 0:
    snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]
    bin_ends = [2, 2.5, 3, 3, 3, 4, 4, 4.5, 4.5, 5, 5]
else:
    snaps = None
    bin_ends = None
print(f'scatter data: {snaps} {bin_ends}')

# Scatter data to all MPI processes
s = comm.scatter(snaps, root=0)
b = comm.scatter(bin_ends, root=0)

# Apply the computation
comm.Barrier()
os.system(f'python3 code/DMhalo/bootstrap.py --snapnum {s} --bin_end {b}')