def add_accret(sim_type, snapnum, snap_idx_table, snap_accret_table):
    import os
    import numpy as np
    from tqdm import tqdm
    import illustris_python as il
    
    # Source data path
    if sim_type == 'DM':
        basePath =  '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L%dn%dTNG'%(205,1250)+'_DM/output/'
    else:
        basePath =  '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L%dn%dTNG'%(205,1250)+'/output/'
    snap_all_mass = il.groupcat.loadHalos(basePath, snapnum, fields='Group_M_Mean200')
    
    # Load halo density profiles
    halos_dir = f'result/DMhalo_density_profiles_phys/TNG300/sim_205_1250_{sim_type}/snap_{snapnum}/final_densities/'
    halos_list = os.listdir(halos_dir)
    # sort the list
    halos_list = sorted(halos_list)
    print(halos_list)
    # Get the bin starts and ends
    bin_starts = [float(halo.split('-')[1])/10 for halo in halos_list]
    bin_ends = [float(halo.split('-')[2].split('.')[0])/10 for halo in halos_list]
    print(bin_starts)
    print(bin_ends)    

    # Load the data
    halo_M, halo_R, bins, densities, accretions = [], [], [], [], []
    for fname, start, end in zip(halos_list, bin_starts, bin_ends):
        print(fname)
        
        data = np.load(os.path.join(halos_dir, fname), allow_pickle=True).item()
        halo_M.append(data['halo_M_Mean200'])
        halo_R.append(data['halo_R_Mean200'])
        bins.append(data['radial_bins'])
        densities.append(data['densities'])
        
        # Get the halo snap indices
        subset_idx = np.where((snap_all_mass >= 10**start) & (snap_all_mass < 10**end))[0]
        print(subset_idx.shape, data['halo_M_Mean200'].shape)
        # Check if the data shape matches
        if subset_idx.shape != data['halo_M_Mean200'].shape:
            # Stop the script
            os._exit(0)
        else:
            pass
        
        # Add the corresponding accretion rates
        add_accret = []
        for idx in tqdm(subset_idx):
            find_idx = np.where(snap_idx_table == idx)[0]
            if len(find_idx) == 0:
                add_accret.append([np.nan])
            elif len(find_idx) == 1:
                add_accret.append([snap_accret_table.iloc[find_idx].values[0]])
            else:
                # Check if the column name matches
                # print(snap_idx_table[find_idx])
                # print(snap_accret_table[find_idx])
                add_accret.append([np.mean(snap_accret_table.iloc[find_idx])])
        add_accret = np.concatenate(add_accret)
        # print(add_accret)
        
        # Accretions
        accretions.append(add_accret)
    
    # Concatenate data
    halo_M = np.concatenate(halo_M)
    halo_R = np.concatenate(halo_R)
    bins = np.concatenate(bins)
    densities = np.concatenate(densities)
    accretions = np.concatenate(accretions)
    print(halo_M.shape, halo_R.shape, bins.shape, densities.shape, accretions.shape)
        
    # Saved dict
    save_dict = {'h': data['h'], 
                'scale_factor': data['scale_factor'], 
                'z': data['z'], 
                'rho_c': data['rho_c'],
                'halo_M_Mean200': halo_M,
                'halo_R_Mean200': halo_R,
                'radial_bins':    bins,
                'densities':      densities,
                'accretion_rate': accretions}
    
    # Create the save path
    save_dir = f'result/DMhalo_density_profiles/TNG300/sim_205_1250_{sim_type}/snap_{snapnum}/final_densities/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    np.save(os.path.join(save_dir, f'bin-{int(bin_starts[0]*10)}-{int(bin_ends[-1]*10)}'), save_dict)
    
    return save_dir

def add_formation_time(sim, snapnum, idx_dict, mass_dict):
    import os
    import h5py
    import numpy as np
    from tqdm import tqdm
    import illustris_python as il
    
    # Load source data
    # -------------------------------------------------------------------------------------------------
    if 'DM' in sim:
        basePath =  '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L%dn%dTNG'%(205,1250)+'_DM/output/'
    elif 'Hydro' in sim:
        basePath =  '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L%dn%dTNG'%(205,1250)+'/output/'
    snap_all_mass = il.groupcat.loadHalos(basePath, snapnum, fields='Group_M_Mean200')
    
    # Load all redshifts
    snap_list = idx_dict.keys()
    z_dict = {}
    for snap in snap_list:
        with h5py.File(il.snapshot.snapPath(basePath, int(snap[5:])), 'r') as f:
            header = dict(f['Header'].attrs.items())
            scale_factor = header['Time']
            z = 1 / scale_factor - 1
            z_dict[snap] = z
    
    # Load halo density profile data
    # -------------------------------------------------------------------------------------------------
    halos_dir = f'result/DMhalo_density_profiles/{sim}/snap_{snapnum}/final_densities/'
    halos_fname = os.listdir(halos_dir)[0]
    data = np.load(os.path.join(halos_dir, halos_fname), allow_pickle=True).item()
    
    # Check halo numbers
    # -------------------------------------------------------------------------------------------------
    # Get the bin starts and ends
    bin_start = float(halos_fname.split('-')[1])/10 
    bin_end = float(halos_fname.split('-')[2].split('.')[0])/10
    print(bin_start, bin_end)
 
    # Get the halo snap indices
    subset_idx = np.where((snap_all_mass >= 10**bin_start) & (snap_all_mass < 10**bin_end))[0]
    print(subset_idx.shape, data['halo_M_Mean200'].shape)
    # Check if the data shape matches
    if subset_idx.shape != data['halo_M_Mean200'].shape:
        # Stop the script
        os._exit(0)
    else:
        pass
    
    # Add the corresponding formation time
    # -------------------------------------------------------------------------------------------------
    # These two snapshots are two early and the results would be meaningless
    if snapnum == 8 or snapnum == 13:
        add_formation_time = np.nan * np.ones(data['halo_M_Mean200'].shape)
    else:
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
                    # Compute the formation time
                    select_mass = {key: value[find_idx] for key, value in mass_dict.items()}
                    match_zs = compute_formation_time(select_mass, snapnum, z_dict)
                    mean_z = np.mean(match_zs)
                    add_formation_time.append([mean_z])
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
    match_indices = np.argmin(abs_diffs, axis=0) 
    match_snaps = [list(mass.keys())[idx] for idx in match_indices]
    match_zs = [z_dict[snap] for snap in match_snaps]
    return match_zs