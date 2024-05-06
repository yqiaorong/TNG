def compt_density_profile(coordinates, haloCM, halo_R_Mean200, 
                          radius_range=[0.01, 5], number_of_bins=85):
    import numpy as np
    radial_coordinates = np.linalg.norm(coordinates - haloCM, axis=1) 
    
    # Create bins
    radius_bins = np.logspace(np.log10(radius_range[0]*halo_R_Mean200), 
                                 np.log10(radius_range[1]*halo_R_Mean200), number_of_bins)
    density_bins = np.empty((number_of_bins))
    
    # Compute the density profile 
    for i, r in enumerate(radius_bins):
        # mass enclosed in bins
        if i == 0:
            mass_bin = len(radial_coordinates[(radial_coordinates <= r)]) 
            volume_bin =  4/3 * np.pi * r**3 
        else:
            mass_bin = len(radial_coordinates[(radius_bins[i-1] < radial_coordinates) & 
                                                    (radial_coordinates <= r)]) 
            volume_bin =  4/3 * np.pi * (r**3 - radius_bins[i-1]**3) 
        # Compute density
        density_bins[i] = mass_bin / volume_bin
    
    return density_bins, radius_bins

def stacked_density_profile(file_list, mass_criteria, DMsoften, use_bootstrap=True):
    
    import numpy as np
    
    # Initialize the output
    scaled_radius_profile, density_profile = [], []
    R200_list = []
    
    ### Load all density profile data ###
    # Iterate over halos
    for file in file_list:
        # load data
        data = np.load(file, allow_pickle=True).item()
        halo_M_Mean200 = data['halo_M_Mean200']
        # apply the mass criteria 
        if (halo_M_Mean200 > 10**mass_criteria[0]) & (halo_M_Mean200 < 10**mass_criteria[1]):
            # radius profile
            halo_R_Mean200 = data['halo_R_Mean200']
            scaled_radius_profile = data['radius_profile']/halo_R_Mean200
            R200_list.append(halo_R_Mean200)
            # density profile
            density_profile.append(data['density_profile'])
        else:
            pass 
    # Get the number of DM halos
    num_halo = np.array(density_profile).shape[0]
    if num_halo == 0:
        print('no halos')
    else:
        print(f'num of halos: {num_halo}')
        
    density_profile = np.array(density_profile).T # shape: (N radii, num_of_halos)
    
    ### Find the median density profile ###
    if use_bootstrap == True:
        medians, errors = [], []
        for rho_data_at_r in density_profile:
            resampled_median_rho_data_at_r = bootstrap(rho_data_at_r, np.median)
            median_with_error = np.percentile(resampled_median_rho_data_at_r, [16, 50, 84])
            # Append new data to the lists
            medians.append(median_with_error[1])
            errors.append([median_with_error[0], median_with_error[2]])
        # Convert the list to array
        medians = np.array(medians) # shape: (N radii, 1)
        errors = np.array(errors) # shape: (N radii, 2)
    else:
        medians = np.percentile(density_profile, 50, axis=1) # shape: (N radii, 1)
        errors = np.percentile(density_profile, [16,84], axis=1).T # shape: (N radii, 2)
    
    ### Compute the minimum radius in percentage in the stack
    if len(R200_list) != 0:
        min_ratio = DMsoften / min(R200_list)
        new_radius_profile = scaled_radius_profile - min_ratio
        min_radius_index = new_radius_profile.tolist().index(min(abs(new_radius_profile)))
    else:
        min_radius_index = 0
    del R200_list
       
    return np.array(scaled_radius_profile[min_radius_index:]), medians[min_radius_index:], errors[min_radius_index:], num_halo

def DM(snapNum):
    
    import h5py
    import illustris_python as il
    
    basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L205n1250TNG/output'
    
    # Groupcat
    Header = il.groupcat.loadHeader(basePath, snapNum)
    Boxsize = Header['BoxSize'] # ckpc/h  
    
    # snapshot
    with h5py.File(il.snapshot.snapPath(basePath, snapNum), 'r') as f:
        header = dict(f['Header'].attrs.items())
        nPart = il.snapshot.getNumPart(header)
        DM_nPart = nPart[1]
        # DM mass
        DMmass = header['MassTable'][1]
    
    soft_length = Boxsize / (DM_nPart**(1/3) * 40) # ckpc/h

    return soft_length, DMmass

def bootstrap(x, statfunc, Nboots=32):
    import numpy as np
    
    x = np.array(x)
    
    resampled_stat = []
    for k in range(Nboots):
        index = np.random.randint(0, len(x), len(x))
        sample = x[index]
        stat_value = statfunc(sample)
        resampled_stat.append(stat_value)
    
    return resampled_stat

def gradient(r, rho):
    import numpy as np
    
    slopes = []
    for i in (range(r.shape[0])):
        if i >= 4:
            slope = (1/12 * np.log(rho[i-4]) - 2/3 * np.log(rho[i-3]) + 
                2/3 * np.log(rho[i-1]) - 1/12 * np.log(rho[i])) / (
                    np.log(r[i]) - np.log(r[i-4]))
            slopes.append(slope)
    
    return np.array(slopes)