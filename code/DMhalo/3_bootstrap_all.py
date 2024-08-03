import os
from mpi4py import MPI

# Parallel computing
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(f'size: {size}, rank: {rank}')

# Define initial data (the data you want to split) on the root MPI process
if rank == 0:
    data = [8, 13, 17, 21, 25, 
            # 33, 40, 50, 67, 78, 99
            ]
else:
    data = None
print(f'scatter data: {data}')

# Scatter data to all MPI processes
s = comm.scatter(data, root=0)

# Apply the computation
comm.Barrier()
os.system(f'python3 code/DMhalo/bootstrap.py --snapnum {s}')