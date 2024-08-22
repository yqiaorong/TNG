import os
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',          default=None,                     type=str)
parser.add_argument('--snapnum',      default=264,                      type=int)
parser.add_argument('--bin_start',    default=1,                        type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',      default=5.5,                      type=float) # [10^{10+x} Msun/h]
parser.add_argument('--save_root_dir',default='DMhalo_density_profiles',type=str)
args = parser.parse_args()

# Set up directory
basePath = f'/virgotng/mpa/MTNG/{args.sim}'
snapnum = args.snapnum

load_dir = os.path.join(basePath, f'snapdir_{snapnum:03d}')
load_list = os.listdir(load_dir)
load_list = [fname for fname in load_list if fname.startswith(f'snapshot_{snapnum:03d}') and fname.endswith('hdf5')]

chunks = range(len(load_list))
for chunk in chunks:
    if not os.path.exists(f'result/{args.save_root_dir}/{args.sim}/snap_{snapnum}/intermediate_densities/'+
                          f'bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}/chunk-{chunk}.npy'):
        os.system('python3 code/DMhalo_tgt/subset_tgt_chunk2.py'+
                f' --sim {args.sim} --snapnum {args.snapnum} --chunk_idx {chunk}'+
                f' --bin_start {args.bin_start} --bin_end {args.bin_end}'
                f' --save_root_dir {args.save_root_dir}')
    else:
        print(f'chunk-{chunk}_bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}.npy: exist')
print(f'All chunk files in snap {args.snapnum} are finished.')