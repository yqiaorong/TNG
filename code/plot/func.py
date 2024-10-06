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

def load_data(dir, snaps, args):
    import numpy as np
    for isnap, snap in enumerate(snaps): # from low z to high z (present)
            
        # Load data
        data = np.load(dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        z = data['z']
        if args.z_or_a == 'a':
            z = 1/ (z+1) # Here z actually means a
        
        mass_bins = data['mass_bins']
        mass_bins = [10**(10+m) for m in mass_bins]
        num_bins = len(mass_bins)
        
        data = data['final_results']
        if 'MTNG' in dir:
            factors = [1000, 1, 1, 1000]
            data[:, args.feat_idx, :] = data[:, args.feat_idx, :]*factors[args.feat_idx] # convert Rsp in [Mpc] to [kpc]
        
        # Concatenate data
        if isnap == 0:
            tot_z = [np.round(z, 3)]*num_bins
            tot_mass_cuts = mass_bins
            tot_data = data
        else:
            tot_z += [np.round(z, 3)]*num_bins
            tot_mass_cuts += mass_bins
            tot_data = np.concatenate((tot_data, data), axis=0)
    
    return np.array(tot_z), np.array(tot_mass_cuts), tot_data

def float_to_str(bin):
    if str(int(bin*10)).endswith('0'):
        return str(int(bin))
    else:
        return str(bin).replace('.', '-')
    
def load_z(dir, snaps):
    import numpy as np
    data = np.load(dir+f'/snap_{min(snaps)}_Rsp_stats.npy', allow_pickle=True).item()
    z_i = np.round(data['z'], 3)
    data = np.load(dir+f'/snap_{max(snaps)}_Rsp_stats.npy', allow_pickle=True).item()
    z_f = np.round(data['z'], 3)
    return z_i, z_f