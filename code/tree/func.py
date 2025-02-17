# def Lifeline2(treedir, chunk_idx):
#     import h5py
#     from tqdm import tqdm
#     import numpy as np
#     import pandas as pd
#     pd.set_option('future.no_silent_downcasting', True)
#     # from itertools import count
#     # iterator = count(0, 1) 
#     from tqdm import tqdm
    
#     with h5py.File(f'{treedir}/trees.{chunk_idx}.hdf5', 'r') as f:
#         # print(f['TreeHalos'].keys())
        
#         # tot_TreeID = np.unique(f['TreeTable/TreeID'][:])      # (Number of trees in this chunk file,)
#         # print(f'The number of Trees in this chunk file {len(tot_TreeID)}')
#         # print('')
#         # GroupNr         = f['TreeHalos/GroupNr']
#         Group_M_Crit200 = f['TreeHalos/Group_M_Crit200'][:]
#         SnapNum         = f['TreeHalos/SnapNum'][:]           # Convert from Dataset to array    
#         # Descendant      = f['TreeHalos/TreeDescendant']       # It gives the index in this chunk file
#         FirstProgenitor = f[f'TreeHalos/TreeFirstProgenitor'][:] # It gives the index in this chunk file
#         # TreeID          = f[f'TreeHalos/TreeID']              # The unique ID of tree
#         SubhaloID       = f[f'TreeHalos/TreeIndex'][:]           # The "unique" ID of subhalo throughout
        
#         subhalo_dict = {f'{snap}_{id}': idx_in_file 
#                         for snap, (idx_in_file, id) in zip(SnapNum, enumerate(SubhaloID))}
        
#         # Select subhalos from the last snap
#         indices_in_file = np.where(SnapNum == 264)[0]  # Unique

#         if indices_in_file.shape[0] == 0:
#             df = pd.DataFrame(index = [f'snap_{x}' for x in range(264, 15, -1)])
#             pass
#         else:
#             # Create the dataframe
#             last_snap_IDs = SubhaloID[indices_in_file]      # The unique subhalo IDs (each corresponds to different snapshots)
#             last_snap_IDs = np.unique(last_snap_IDs)        # The unique subhalo IDs
            
#             df = pd.DataFrame(index = [f'snap_{x}' for x in range(264, 15, -1)],
#                               columns = last_snap_IDs)
            
#             # Update the last snap
#             FP_IDs = last_snap_IDs
#             del last_snap_IDs
#             for snap in tqdm(range(264, 15, -1)):

#                 FP_indices = np.array([subhalo_dict.get(f'{snap}_{id}', -1) for id in FP_IDs])
#                 df.loc[f'snap_{snap}'] = FP_indices
#                 # Update variables, moving to one prev snap
#                 FP_IDs = FirstProgenitor[FP_indices]

#             # Modify incorrect FP ids
#             df = df.apply(lambda col: np.where(col.cummin() == -1, -1, col))
            
#             # Assign mass to FP_indices
#             df = df.map(lambda x: Group_M_Crit200[x] if x != -1 else np.nan) 

#     return df

# def Lifeline(treedir, chunk_idx):
#     import h5py
#     from tqdm import tqdm
#     import numpy as np
#     import pandas as pd
#     pd.set_option('future.no_silent_downcasting', True)
#     # from itertools import count
#     # iterator = count(0, 1) 
#     from tqdm import tqdm
    
#     with h5py.File(f'{treedir}/trees.{chunk_idx}.hdf5', 'r') as f:
#         # print(f['TreeHalos'].keys())
        
#         # tot_TreeID = np.unique(f['TreeTable/TreeID'][:])      # (Number of trees in this chunk file,)
#         # print(f'The number of Trees in this chunk file {len(tot_TreeID)}')
#         # print('')
#         # GroupNr         = f['TreeHalos/GroupNr']
#         Group_M_Crit200 = f['TreeHalos/Group_M_Crit200'][:]
#         SnapNum         = f['TreeHalos/SnapNum'][:]           # Convert from Dataset to array    
#         # Descendant      = f['TreeHalos/TreeDescendant']       # It gives the index in this chunk file
#         FirstProgenitor = f[f'TreeHalos/TreeFirstProgenitor'][:] # It gives the index in this chunk file
#         # TreeID          = f[f'TreeHalos/TreeID']              # The unique ID of tree
#         SubhaloID       = f[f'TreeHalos/TreeIndex'][:]           # The "unique" ID of subhalo throughout
#         subhalo_dict = {id: idx for idx, id in enumerate(SubhaloID)}
#         subhalo_dict[-1] = -1
        
#         # Select subhalos from the last snap
#         last_snap = 264
        
#         indices_in_file = np.where(SnapNum == last_snap)[0]  # Unique

#         if indices_in_file.shape[0] == 0:
#             df = pd.DataFrame(index = [f'snap_{x}' for x in range(264, 15, -1)])
#             pass
#         else:
#             # Create the dataframe
#             last_snap_IDs = SubhaloID[indices_in_file]      # The unique subhalo IDs
#             last_snap_IDs = np.unique(last_snap_IDs)

#             df = pd.DataFrame(index = [f'snap_{x}' for x in range(264, 15, -1)],
#                               columns = last_snap_IDs)
            
#             # Update the last snap
#             snap = last_snap
#             df.loc[f'snap_{snap}', last_snap_IDs] = last_snap_IDs # unique so far

#             # Update variables, moving to one prev snap
#             FP_indices = np.array([subhalo_dict[id] for id in last_snap_IDs]) # The shape matches now

#             for snap in tqdm(range(263, 15, -1)):
                
#                 FP_IDs = FirstProgenitor[FP_indices]
#                 # print(FP_IDs.shape, FP_indices.shape)
                
#                 # print(df.loc[f'snap_{snap}'].shape, FP_IDs.shape)
#                 df.loc[f'snap_{snap}'] = FP_IDs

#                 FP_indices = np.array([subhalo_dict.get(id, -1) for id in FP_IDs])

#             # Modify incorrect FP ids
#             df = df.apply(lambda col: np.where(col.cummin() == -1, -1, col))
            
#             # Assign mass to IDs
#             subid_to_mass = dict(zip(SubhaloID, Group_M_Crit200))
#             subid_to_mass[-1] = np.nan
#             df = df.map(lambda x: subid_to_mass.get(x, x))

#     return df

def calc_tdyn(basePath, snap_dict, current_snap):
    import h5py
    import numpy as np
    import astropy.units as u
    import illustris_python as il
    from astropy.cosmology import Planck15, z_at_value, FlatLambdaCDM
    
    # # Load redshifts in the snap_list
    # z_list = []
    # for snap in snap_list:
    #     with h5py.File(il.snapshot.snapPath(basePath, int(snap[5:])), 'r') as f:
    #         header = dict(f['Header'].attrs.items())
    #         scale_factor = header['Time']
    #         z = 1 / scale_factor - 1
    #     z_list.append(z)
    
    # Load the current snapshot redshift
    with h5py.File(il.snapshot.snapPath(basePath, current_snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        current_z = 1 / scale_factor - 1

    # Define own cosmology
    h = 0.6774
    Om0 = 0.3089
    Ob0 = 0.0486
    cosmo = FlatLambdaCDM(H0 = h * 100 * u.km / u.s / u.Mpc, 
                            Om0=Om0, Ob0=Ob0, Tcmb0=2.725)
    
    t = cosmo.age(current_z) # The cosmological time at current snapshot / z
    H = cosmo.H(current_z)   # The Hubble parameter at current snapshot / z
    
    # Cosmological dynamic time
    t_H   = 1 / H
    t_dyn = t_H / (5 * np.sqrt(Om0))
    t_dyn = t_dyn.to('yr').value * 1E-9 # t_dyn in unit Gyr
    print(f"One dynamical time: {t_dyn} Gyr at current z = {current_z}")
    
    # Find z one dynamical time before current z
    p_z = z_at_value(cosmo.age,  t - t_dyn * u.Gyr) # t_dyn in unit Gyr
    prev_snap = min(snap_dict, key=lambda snap: abs(snap_dict[snap] - p_z))
    prev_z = snap_dict[prev_snap]
    
    # Find z one dynamical time after current z
    if np.round(current_z, 2) < 0.5:
        l_z = 0.00
    else:
        l_z = z_at_value(cosmo.age, t + t_dyn * u.Gyr) # t_dyn in unit Gyr
    later_snap = min(snap_dict, key=lambda snap: abs(snap_dict[snap] - l_z))
    later_z = snap_dict[later_snap]

    print('The previous redshift: ', prev_snap, prev_z)
    print('The later redshift: ', later_snap, later_z)
    return prev_snap, later_snap