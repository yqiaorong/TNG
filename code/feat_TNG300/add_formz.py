import os
import numpy as np
import argparse
from tqdm import tqdm
import illustris_python as il
from func_accret import *
from colossus.cosmology import cosmology
from colossus.lss import peaks
cosmology.setCosmology('planck15')
from mpi4py import MPI

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default='TNG300/sim_205_1250_Hydro', type=str)
parser.add_argument('--snapnum', default=None, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add formation time to halo data <<<')
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

    # Find formation time from peak height
    formz = []
    for current_mass, subhalo_id in tqdm(
                                        zip(data['halo_M_Mean200'][start:end], data['FirstSub'][start:end]),
                                        total = end-start, desc = f"Rank {rank}"
                                        ): # Physical mass!
        
        # Find the first subhalo history
        result = il.lhalotree.loadTree(basePath, args.snapnum, subhalo_id, fields=['SnapNum','SubhaloGrNr'], onlyMPB=True)
        if result is None:
            formz.append([np.nan])
        else:
            # Find peak height one by one
            result['pH'] = []
            for s, gr_id in zip(result['SnapNum'], result['SubhaloGrNr']):
                halo_info = il.groupcat.loadSingle(basePath, s, haloID=gr_id)
                halo_info = halo_info['Group_M_Mean200'] # Comoving mass!
                if halo_info == 0:
                    result['pH'].append(np.nan)
                else:
                    ph = peaks.peakHeight(halo_info*10**10/data['h'], z_dict[s]) # physical mass!
                    result['pH'].append(ph)
            
            # Find the index in result that has peak_height closest to 1 
            mask = (np.array(result['pH']) >= 0.9) & (np.array(result['pH']) <= 1.1)
            if not np.any(mask):
                formz.append([np.nan])
            else:
                closest_idx = np.argmin(np.abs(np.array(result['pH'])[mask] - 1.0))
                final_idx = np.where(mask)[0][closest_idx]
                form_snap = result['SnapNum'][final_idx]
                form_z = z_dict[form_snap]
                # print(form_snap, form_z, result['pH'][final_idx])
                formz.append([form_z])
    formz_local = np.array(formz)
    
    # Gather results on the root process
    formz_all = comm.gather(formz_local, root=0)

    if rank == 0:
        formz_full = np.concatenate(formz_all, axis=0)  
        data['formz'] = formz_full.squeeze()
        print('Final formz shape:', formz_full.squeeze().shape)
        print(data.keys())
        np.save(halos_dir + fname, data)
        print('Saved!')
             
    # print(np.concatenate(formz).shape)
    # # Add accretion rates to data
    # data['formz'] = np.concatenate(formz)
    # print(data.keys())
    # # Save data
    # np.save(halos_dir+fname, data)
    # print('Saved! ')