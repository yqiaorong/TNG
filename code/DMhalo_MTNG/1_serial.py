import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sim',             default=None, type=str)
parser.add_argument('--snapnum',         default=None, type=int)
parser.add_argument('--start_chunk_idx', default=None, type=int)
parser.add_argument('--bin_start',       default=None, type=float)
parser.add_argument('--bin_end',         default=None, type=float)
args = parser.parse_args()

chunk_size = 4
start = args.start_chunk_idx
parent_dir = f'/nfs/mvogelsblab001/Users/s_qyu/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/'

for i in range(start, int(start+chunk_size)):
    if not os.path.exists(f'{parent_dir}/intermediate_densities/bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}/chunk-{i}.npy'):
        os.system(f'python3.11 code/DMhalo/subset_tgt_chunk-gpu.py --sim {args.sim} --snapnum {args.snapnum} '+
                  f'--bin_start {args.bin_start} --bin_end {args.bin_end} --chunk_idx {i}')