def Lifeline(treedir, chunk_idx):
    import h5py
    from tqdm import tqdm
    import numpy as np
    import pandas as pd
    pd.set_option('future.no_silent_downcasting', True)
    # from itertools import count
    # iterator = count(0, 1) 
    from tqdm import tqdm
    
    with h5py.File(f'{treedir}/trees.{chunk_idx}.hdf5', 'r') as f:
        # print(f['TreeHalos'].keys())
        
        # tot_TreeID = np.unique(f['TreeTable/TreeID'][:])      # (Number of trees in this chunk file,)
        # print(f'The number of Trees in this chunk file {len(tot_TreeID)}')
        # print('')
        # GroupNr         = f['TreeHalos/GroupNr']
        Group_M_Crit200 = f['TreeHalos/Group_M_Crit200'][:]
        SnapNum         = f['TreeHalos/SnapNum'][:]           # Convert from Dataset to array    
        # Descendant      = f['TreeHalos/TreeDescendant']       # It gives the index in this chunk file
        FirstProgenitor = f[f'TreeHalos/TreeFirstProgenitor'][:] # It gives the index in this chunk file
        # TreeID          = f[f'TreeHalos/TreeID']              # The unique ID of tree
        SubhaloID       = f[f'TreeHalos/TreeIndex'][:]           # The "unique" ID of subhalo throughout
        subhalo_dict = {id: idx for idx, id in enumerate(SubhaloID)}
        subhalo_dict[-1] = -1
        
        # Select subhalos from the last snap
        last_snap = 264
        
        indices_in_file = np.where(SnapNum == last_snap)[0]  # Unique

        if indices_in_file.shape[0] == 0:
            df = pd.DataFrame(index = [f'snap_{x}' for x in range(264, 15, -1)])
            pass
        else:
            # Create the dataframe
            last_snap_IDs = SubhaloID[indices_in_file]      # The unique subhalo IDs
            last_snap_IDs = np.unique(last_snap_IDs)

            df = pd.DataFrame(index = [f'snap_{x}' for x in range(264, 15, -1)],
                              columns = last_snap_IDs)
            
            # Update the last snap
            snap = last_snap
            df.loc[f'snap_{snap}', last_snap_IDs] = last_snap_IDs # unique so far

            # Update variables, moving to one prev snap
            FP_indices = np.array([subhalo_dict[id] for id in last_snap_IDs]) # The shape matches now

            for snap in tqdm(range(263, 15, -1)):
                
                FP_IDs = FirstProgenitor[FP_indices]
                # print(FP_IDs.shape, FP_indices.shape)
                
                # print(df.loc[f'snap_{snap}'].shape, FP_IDs.shape)
                df.loc[f'snap_{snap}'] = FP_IDs

                FP_indices = np.array([subhalo_dict.get(id, -1) for id in FP_IDs])

            # Modify incorrect FP ids
            df = df.apply(lambda col: np.where(col.cummin() == -1, -1, col))
            
            # Assign mass to IDs
            subid_to_mass = dict(zip(SubhaloID, Group_M_Crit200))
            subid_to_mass[-1] = np.nan
            df = df.map(lambda x: subid_to_mass.get(x, x))

    return df

# def lifeline2(treedir, chunk_idx):
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
        
#         tot_TreeID = np.unique(f['TreeTable/TreeID'][:])      # (Number of trees in this chunk file,)
#         print(f'The number of Trees in this chunk file {len(tot_TreeID)}')
#         print('')
#         # GroupNr         = f['TreeHalos/GroupNr']
#         Group_M_Crit200 = f['TreeHalos/Group_M_Crit200'][:]
#         SnapNum         = f['TreeHalos/SnapNum'][:]           # Convert from Dataset to array    
#         # Descendant      = f['TreeHalos/TreeDescendant']       # It gives the index in this chunk file
#         FirstProgenitor = f[f'TreeHalos/TreeFirstProgenitor'][:] # It gives the index in this chunk file
#         TreeID          = f[f'TreeHalos/TreeID']              # The unique ID of tree
#         SubhaloID       = f[f'TreeHalos/TreeIndex'][:]          # The unique ID of subhalo throughout
        
#         # Select subhalos from the last snap
#         last_snap = 264
        

#         indices_in_tree = np.where(SnapNum == last_snap)[0]  
#         print(indices_in_tree.shape)
#         # Sanitiy check of treeID
#         # print(f'Sanity check: selected {indices_in_tree.shape} subhalos belong to tree {np.unique(TrID[indices_in_tree])}')

#         if indices_in_tree.shape[0] == 0:
#             pass
#         else:
#             # Create the dataframe
#             start_IDs = SubhaloID[indices_in_tree]      # The unique subhalo IDs
#             FP_IDs = start_IDs
#             df = np.empty((len(range(16, 265)), start_IDs.shape[0]))
#             # df = pd.DataFrame(index = [f'snap_{x}' for x in range(16, 265)], columns = column_names)
#             print(df.shape)
#             # Update the last snap
#             snap = last_snap
#             isnap = 0
#             # print(np.all(df.columns == column_names))
#             # print(FP_IDs.shape)
#             # print(df.loc[f'snap_264', df.columns].shape)
#             # df.loc[f'snap_264', column_names] = FP_IDs
            
            
            
#             # Update variables, moving to one prev snap
#             FP_indices = indices_in_tree
#             df[isnap] = FP_IDs
#             print('FP indices', FP_indices.shape)
            
#             assign_col_idx = np.where(FP_IDs = start_IDs)
            
#             for snap in tqdm(range(263, 15, -1)):
                
#                 FP_IDs = FirstProgenitor[FP_indices]
#                 # print('FP IDs', FP_IDs.shape, np.unique(FP_IDs).shape, min(FP_IDs))
#                 # print(assign_column_names.shape, FP_IDs.shape)
#                 if assign_column_names.shape == FP_IDs.shape:
#                     print('true')
#                     print(df.loc[f'snap_{snap}', assign_column_names].shape)
#                     df.loc[f'snap_{snap}', assign_column_names] = np.array(FP_IDs)
#                 else:
#                     print('pass')

#                 print('column names', column_names.shape, 'FP IDs', FP_IDs.shape)
        
#                 end_columns = df.columns[df.loc[f'snap_{snap}'] == -1]
#                 print(np.unique(column_names.shape), np.unique(end_columns.shape))
#                 df.loc[:, end_columns] = df.loc[:, end_columns].fillna(-1)# .infer_objects(copy=False)
#                 print('end columns', end_columns.shape)
                
#                 # print(np.all(np.isin(column_names, end_columns)))
#                 # column_names = np.array([col for col in column_names if col not in end_columns])
#                 assign_column_names = np.setdiff1d(column_names, end_columns)
#                 print('assign column names', assign_column_names.shape)
#                 FP_indices = np.where(np.isin(SubhaloID, FP_IDs))[0]
#                 print('FP indices', FP_indices.shape)
#                 print('')
                
#             # Assign mass to IDs
#             subid_to_mass = dict(zip(SubhaloID, Group_M_Crit200))
#             df = df.map(lambda x: subid_to_mass.get(x, x))
#     return df

# def get_field_values_of_lifeline(basePath, df, group_field=None, coords_idx=None):
#     """Only specify the coordinates idx of group_fields includes GroupPos"""
#     from tqdm import tqdm
#     import numpy as np
#     import pandas as pd
#     import illustris_python as il
    
#     for irow, row in tqdm(df.iterrows(), desc='snaps'):

#         snapnum = int(irow[5:])
#         print(irow)
#         # Load halo field
#         Halos = il.groupcat.loadHalos(basePath, snapnum, fields=group_field)
#         if group_field == 'GroupPos':
#             Halos = Halos[:, coords_idx]
#         else:
#             pass
        
#         # Load subhalo global index 
#         # SubhaloGrNr = il.groupcat.loadSubhalos(basePath, snapnum, fields='SubhaloGroupNr') # It gives the "index" into the groups catalogue
        
#         # Substitute the subhalo global index with the parent group_field
#         df.loc[irow] = [Halos[int(idx)] if pd.notna(idx) else np.nan for idx in row]
    
#     return df  