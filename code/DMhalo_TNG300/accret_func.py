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