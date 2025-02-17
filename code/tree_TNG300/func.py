# This function is quoted from Thesan offsets documentation so never change any of it.

def get_global_index(simpath, snap, offset_type, local_index, chunk):

    """
    Returns the global index (i.e. across all chunks, considered sequentially) of a particle/group/subgroup.

    Parameters
    ----------
    simpath      : string
                   path to the base simulation directory
    snap         : int
                   number of the snapshot to consider
    offset_type  : string
                   string identifying the type of object considered. Accepted values: 'particle', 'group', 'subhalo'
    local_index  : int
                   index of the object in the chunk file
    chunk        : int
                   number of the chunk file the object resides in

    Returns
    -------
    global_index  : int
                    index of the object across all chunk files
    """
    
    import h5py

    assert(offset_type in ['particle', 'group', 'subhalo'])
    if offset_type == 'particle':
        offset_file_key = 'FileOffsets/SnapByType'
    elif offset_type == 'group':
        offset_file_key = 'FileOffsets/Group'
    elif offset_type == 'subhalo':
        offset_file_key = 'FileOffsets/Subhalo'

    with h5py.File(f'{simpath}/postprocessing/offsets/offsets_{snap:03d}.hdf5', 'r') as offset_file:
        global_index = offset_file[offset_file_key][chunk] + local_index

    return global_index


def get_chunk_and_local_index(simpath, snap, offset_type, global_index, ptype=-1):

    """
    Returns the chunk file number and local index (i.e. within the chunk file) of a particle/group/subgroup.

    Parameters
    ----------
    simpath       : string
                    path to the base simulation directory
    snap          : int
                    number of the snapshot to consider
    offset_type   : string
                    string identifying the type of object considered. Accepted values: 'particle', 'group', 'subhalo'
    global_index  : int
                    index of the object across all chunk files
    ptype         : int, required for offset_type == 'particle', ignored otherwise
                  : particle type

    Returns
    -------
    chunk        : int
                   number of the first chunk file the object resides in
    local_index  : int
                   index of the object in the chunk file
    """
    
    import h5py
    import numpy as np
    
    assert(offset_type in ['particle', 'group', 'subhalo'])
    if offset_type == 'particle':
        assert(ptype>=0 and ptype<6)
    
    
    if offset_type == 'particle':
        offset_file_key = 'FileOffsets/SnapByType'
    elif offset_type == 'group':
        offset_file_key = 'FileOffsets/Group'
    elif offset_type == 'subhalo':
        offset_file_key = 'FileOffsets/Subhalo'

    with h5py.File(f'{simpath}/postprocessing/offsets/offsets_{snap:03d}.hdf5', 'r') as offset_file:
        if offset_type == 'particle':
            chunk = np.where(offset_file[offset_file_key][:, ptype] <= global_index)[0][-1]
            local_index = global_index - offset_file[offset_file_key][chunk, ptype]
        else:
            chunk = np.where(offset_file[offset_file_key][()]       <= global_index)[0][-1]
            local_index = global_index - offset_file[offset_file_key][chunk]
        print(offset_file[offset_file_key])
    return chunk, local_index


def lifeline(parent_dir, subhalo_global_idx, df, last_snap_idx=99):
    import os
    import h5py
    import pandas as pd
    from itertools import count
    
    ### Find subhalo with global index's position in MergerTree ###
    # Offsets
    offset_path = os.path.join(parent_dir, f'postprocessing/offsets/offsets_{last_snap_idx:03d}.hdf5')
    with h5py.File(offset_path, 'r') as f:
        
        File = f['Subhalo/LHaloTree/File']
        Num = f['Subhalo/LHaloTree/Num']
        Index = f['Subhalo/LHaloTree/Index']
        
        tree_chunk_idx = File[subhalo_global_idx]
        treeX = Num[subhalo_global_idx]
        subhalo_intreeX_idx = Index[subhalo_global_idx]
        
        # print(f'Subhalo gobal index {subhalo_global_idx} is stored in tree chunk file index: {tree_chunk_idx}.')
        # print(f'In tree chunk file index {tree_chunk_idx}, the subhalo is stored in Tree{treeX}.')
        # print(f'In Tree{treeX}, the subhalo has index: {subhalo_intreeX_idx}.')
        
    if tree_chunk_idx == -1:
        pass
    else: 
        # enter the global index of current subhalo into the dataframe
        df.loc[f'snap_{last_snap_idx}', f'global_idx_{subhalo_global_idx}'] = subhalo_global_idx
        
        ### Finding the lifeline of global index ###
        # Tree 
        tree_path = os.path.join(parent_dir,f'postprocessing/trees/LHaloTree/trees_sf1_{last_snap_idx:03d}.{tree_chunk_idx}.hdf5')
        with h5py.File(tree_path, 'r') as f:
        
            # In treeX
            SubhaloNumber = f[f'Tree{treeX}/SubhaloNumber'] # This is the global index of the subhalo at corresponding snapshot
            SnapNum = f[f'Tree{treeX}/SnapNum']
            Descendant = f[f'Tree{treeX}/Descendant']
            FirstProgenitor = f[f'Tree{treeX}/FirstProgenitor']
            
            # Initials: these are the index of first progenitor and descendant within TreeX
            FP_idx = FirstProgenitor[subhalo_intreeX_idx]
            D_idx = Descendant[subhalo_intreeX_idx]
            
            # Set up the iterator 
            iterator = count(0, 1)
            
            # First progenitor
            for item in iterator:
                if FP_idx == -1:
                    break
                df.loc[f'snap_{SnapNum[FP_idx]}', f'global_idx_{subhalo_global_idx}'] = SubhaloNumber[FP_idx]
                # Update FP_idx
                FP_idx = FirstProgenitor[FP_idx]
                
            # Descendant
            for item in iterator:
                if D_idx == -1:
                    break
                df.loc[f'snap_{SnapNum[D_idx]}', f'global_idx_{subhalo_global_idx}'] = SubhaloNumber[D_idx]
                # Update D_idx
                D_idx = Descendant[D_idx]
    pass
    # return df


def get_field_values_of_lifeline(simpath, df, group_field=None, coords_idx=None):
    """Only specify the coordinates idx of group_fields includes GroupPos"""
    from tqdm import tqdm
    import numpy as np
    import pandas as pd
    import illustris_python as il
    
    for irow, row in tqdm(df.iterrows(), desc='snaps'):

        snapnum = int(irow[5:])
        print(irow)
    
        # Load halo field
        Halos = il.groupcat.loadHalos(simpath+'output', snapnum, fields=group_field)
        if group_field == 'GroupPos':
            Halos = Halos[:, coords_idx]
        else:
            pass
        
        # Load subhalo global index 
        SubhaloGrNr = il.groupcat.loadSubhalos(simpath+'output', snapnum, fields='SubhaloGrNr')

        # Substitute the subhalo global index with the parent halo mass
        df.loc[irow] = [Halos[SubhaloGrNr[int(idx)]] if pd.notna(idx) else np.nan for idx in row]
    
    return df

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