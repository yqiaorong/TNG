import os
import h5py
import numpy as np
import argparse
from tqdm import tqdm
import illustris_python as il
from func_accret import *
from mpi4py import MPI

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default='TNG300/sim_205_1250_Hydro', type=str)
parser.add_argument('--snapnum', default=17, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add mergerz to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Create a snapshot and redshift dict
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if 'DM' in args.sim:
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'_DM/output/'
elif 'Hydro' in args.sim:
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'/output/'
    
z_dict, a_dict = {}, {}
for snap in range(100):
    header = il.groupcat.loadHeader(basePath, snap)
    z_dict[snap] = header['Redshift']
    a_dict[snap] = header['Time']



# Load halo data
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    # Load data
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    print(data.keys())
    print(data['halo_M_Mean200'].shape)
    
    # Parallel computing
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    print(f'size: {size}, rank: {rank}')
    
    if rank == 0:
        split_indices = np.linspace(0, len(data['FirstSub']), size + 1).astype(int)
    else:
        split_indices = None
    split_indices = comm.bcast(split_indices, root=0)

    start = split_indices[rank]
    end = split_indices[rank + 1]
    print(f'Rank {rank} scatter indices: {start} - {end}')

    # Apply the computation
    comm.Barrier()
    
    # Find last merger time
    mergerz = []
    # for current_mass, subhalo_id in tqdm(
    #                                     zip(data['halo_M_Mean200'][start:end], data['FirstSub'][start:end]),
    #                                     total = end-start, desc = f"Rank {rank}" # Physical mass!
    #                                     ):
    for current_mass, subhalo_id in tqdm(zip(data['halo_M_Mean200'][:5], data['FirstSub'][:5]), 
                                         total=len(data['FirstSub'])):
        # Find the first subhalo history
        result = il.lhalotree.loadTree(basePath, args.snapnum, subhalo_id, 
                                       fields=['SnapNum','SubhaloNumber','SubhaloGrNr','NextProgenitor'], 
                                       onlyMPB=True)

        if result is None:
            mergerz.append([np.nan])
        else:
            # Filter out the result with NP = -1
            mask = result['NextProgenitor'] != -1
            f_result = {k: v[mask] if isinstance(v, np.ndarray) and len(v) == len(mask) else v
                            for k, v in result.items()}
            del result
            print(f_result)
            
            # Find tree
            f_result['NPmass'], f_result['Grmass'] = [], []
            for sub_s, sub_gr_id, sub_id, NPIndex in zip(f_result['SnapNum'], f_result['SubhaloGrNr'],
                                                         f_result['SubhaloNumber'], f_result['NextProgenitor']):
                TreeFile, TreeIndex, TreeNum = il.lhalotree.treeOffsets(basePath, sub_s, sub_id)
                print('tree:', TreeFile, TreeIndex, TreeNum)
                if TreeFile == -1:
                    print('No tree found!')
                    f_result['NPmass'].append(np.nan)
                else:
                    # Load full next progenitor data
                    with h5py.File(il.lhalotree.treePath(basePath, TreeFile), 'r') as f:
                        tree = f[f'Tree{TreeNum}'] 
                        
                        snapnums = tree['SnapNum'][:]
                        group_ids = tree['SubhaloGrNr'][:]
                        
                        # Find the subhaloID of the NPIndex
                        NP_s, NP_gr_id = snapnums[NPIndex], group_ids[NPIndex]
                        print(NP_gr_id, NP_s)
                        
                        # Find the NP group mass
                        NP_mass = il.groupcat.loadSingle(basePath, NP_s, haloID=NP_gr_id)
                        NP_mass = NP_mass['Group_M_Mean200'] # Comoving mass!
                        
                        # Find the current group mass
                        sub_mass = il.groupcat.loadSingle(basePath, sub_s, haloID=sub_gr_id)
                        sub_mass = sub_mass['Group_M_Mean200'] # Comoving mass!
                        
                        f_result['NPmass'].append(NP_mass)
                        f_result['Grmass'].append(sub_mass)
            print(f_result)
            print('')
            
            # Find the mas ratio
            f_result['mass_ratio'] = np.array(np.array(f_result['NPmass']) / np.array(f_result['Grmass']))
            # Find the index > 0.1
            valid_idx = np.where(f_result['mass_ratio'] > 0.1)[0]
            print(valid_idx)
            if len(valid_idx) == 0:
                mergerz.append([np.nan])
            else:
                # Find the first idx
                last_idx = valid_idx[0]
                # Find the corresponding snapnum
                merger_snap = f_result['SnapNum'][last_idx]+1
                merger_z = z_dict[merger_snap]
                mergerz.append([merger_z]) # Merge at descendant snapshot
     
    mergerz = np.array(mergerz)
    print(mergerz) 
    print(np.concatenate(mergerz).shape)
    # # Add accretion rates to data
    # data['mergerz'] = np.concatenate(mergerz)
    # # Save data
    # np.save(halos_dir+fname, data)
    # print('Saved! ')