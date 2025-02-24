def add_formation_time(args, idx_dict, mass_dict, submass_dict):
    import os
    import h5py
    import numpy as np
    from tqdm import tqdm
    import illustris_python as il
    
    sim, snapnum = args.sim, args.snapnum
    
    # Load source data
    # -------------------------------------------------------------------------------------------------
    snap_all_mass = il.groupcat.loadHalos(f'/virgotng/mpa/MTNG/{sim}/output/', snapnum, fields='Group_M_Mean200')
    
    # Load all redshifts
    snap_list = idx_dict.keys()
    z_dict = {}
    for snap in snap_list:
        with h5py.File(il.snapshot.snapPath(f'/virgotng/mpa/MTNG/{sim}/output/', int(snap[5:])), 'r') as f:
            header = dict(f['Header'].attrs.items())
            scale_factor = header['Time']
            z = 1 / scale_factor - 1
            z_dict[snap] = z
    
    # Load halo density profile data
    # -------------------------------------------------------------------------------------------------
    halos_dir = f'result/DMhalo_density_profiles/MTNG/{sim}/snap_{snapnum}/final_densities/'
    halos_fname = os.listdir(halos_dir)[0]
    data = np.load(halos_dir+halos_fname, allow_pickle=True).item()
    
    # Check halo numbers
    # -------------------------------------------------------------------------------------------------
    # Get the bin starts and ends
    bin_start = float(halos_fname.split('-')[1])/10 
    bin_end = float(halos_fname.split('-')[2].split('.')[0])/10
    print(bin_start, bin_end)
 
    # Get the halo snap indices
    subset_idx = np.where((snap_all_mass >= 10**bin_start) & (snap_all_mass < 10**bin_end))[0]
    print(subset_idx.shape, data['halo_M_Mean200'].shape)
        
    # Add the corresponding formation time
    # -------------------------------------------------------------------------------------------------
    add_formation_time = []
    for idx in tqdm(subset_idx):
        find_idx = np.where(idx_dict[f'snap_{snapnum}'] == idx)[0]
        
        if len(find_idx) == 0:
            # print('no halo found')
            add_formation_time.append([np.nan])
        else:
            # Check if mass matches
            if np.all(np.abs(mass_dict[f'snap_{snapnum}'][find_idx] - snap_all_mass[idx]) < 0.01):
                # print('mass matches')
                if len(find_idx) == 1:
                    find_idx = find_idx[0]
                else:
                    # Select the one with the most massive subhalo mass
                    find_idx = find_idx[np.argmax(submass_dict[f'snap_{snapnum}'][find_idx])]
  
                # Compute the formation time
                select_mass = {key: value[find_idx] for key, value in mass_dict.items()}
                match_z = compute_formation_time(select_mass, snapnum, z_dict)
                add_formation_time.append([match_z])
            else:
                # print('mass does not match')
                exit()
    add_formation_time = np.concatenate(add_formation_time)
        
    print(add_formation_time.shape, data['halo_M_Mean200'].shape)
    
    # Add formation time to halo data
    # -------------------------------------------------------------------------------------------------
    data['formation_time'] = add_formation_time
    print(data.keys())
    np.save(f'{halos_dir}/{halos_fname}', data)
    print('data saved')
    
def compute_formation_time(mass, snapnum, z_dict):
    
    """
    Find the row index in `column_mass` where the column values are equal to 
    half of the values in row `snap_{snapnum}`.
    
    Parameters:
    mass : dict
        A dictionary where keys are 'snap_{snapnum}' and values are 1D NumPy arrays.
    snapnum : int
        The current snapshot.
    
    Returns:
    dict
        A dictionary where keys are column names and values are lists of row indices
        that satisfy the condition.
    """
    import numpy as np
    mass_array = np.array([mass[key] for key in mass.keys()]) 
    half_mass = mass[f'snap_{snapnum}'] / 2 
    abs_diffs = np.abs(mass_array - half_mass)
    # Find the closest matching row (snapshot) for each object
    match_idx = np.argmin(abs_diffs, axis=0) 
    match_snap = list(mass.keys())[match_idx] 
    match_zs = z_dict[match_snap] 
    return match_zs