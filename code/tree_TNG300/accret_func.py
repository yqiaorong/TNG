def add_accret(sim_type, snapnum, snap_idx_table, snap_mass_table, snap_accret_table, snap_subhalo_mass_table):
    import os
    import numpy as np
    from tqdm import tqdm
    import illustris_python as il
    
    # Source data path
    if sim_type == 'DM':
        basePath =  '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L%dn%dTNG'%(205,1250)+'_DM/output/'
    else:
        basePath =  '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L%dn%dTNG'%(205,1250)+'/output/'
    snap_all_mass = il.groupcat.loadHalos(basePath, snapnum, fields='Group_M_Mean200') # Comoving mass!!!
    
    # Load halo density profiles
    halos_dir = f'result/DMhalo_density_profiles/TNG300/sim_205_1250_{sim_type}/snap_{snapnum}/final_densities/'
    halos_list = os.listdir(halos_dir)
    print(halos_list)
    # Get the bin starts and ends
    bin_starts = [float(halo.split('-')[1])/10 for halo in halos_list]
    bin_ends = [float(halo.split('-')[2].split('.')[0])/10 for halo in halos_list]
    print(bin_starts)
    print(bin_ends)    

    # Load the data
    for fname, start, end in zip(halos_list, bin_starts, bin_ends):
        print(fname)
        
        data = np.load(os.path.join(halos_dir, fname), allow_pickle=True).item()
        
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
            # Double check if the mass matches
            if np.all(np.abs(snap_mass_table.iloc[find_idx].values - snap_all_mass[idx]) < 0.001) == False:
                print('mass does not match')
                exit()
            else:
                if len(find_idx) == 0:
                    # print('no halo found')
                    add_accret.append([np.nan])
                elif len(find_idx) == 1:
                    add_accret.append([snap_accret_table.iloc[find_idx].values[0]])
                else:
                    # Choose the one with most massive subhalo
                    right_idx_in_find_idx = np.where(snap_subhalo_mass_table.iloc[find_idx].values == np.max(snap_subhalo_mass_table.iloc[find_idx].values))[0]
                    # print('idx:', find_idx, right_idx_in_find_idx)
                    # print('subhalo mass:', snap_subhalo_mass_table.iloc[find_idx].values, snap_subhalo_mass_table.iloc[find_idx[right_idx_in_find_idx]].values)
                    add_accret.append([snap_accret_table.iloc[find_idx[right_idx_in_find_idx]].values[0]])
        
        # Accretions
        add_accret = np.array(np.concatenate(add_accret)).squeeze()
        # Print the number of non NaN entries in add_accret
        print(np.sum(~np.isnan(add_accret)), 'are non-NaN in ', add_accret.shape, data['halo_M_Mean200'].shape)
        data['accretions'] = add_accret
        print(data.keys())
        np.save(os.path.join(halos_dir, fname), data)
    
    return 'done'