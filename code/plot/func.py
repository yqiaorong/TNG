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
        # if 'MTNG' in dir:
        #     factors = [1000, 1, 1, 1000]
        #     data[:, args.feat_idx, :] = data[:, args.feat_idx, :]*factors[args.feat_idx] # convert Rsp in [Mpc] to [kpc]
        
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

### Used in plot depth vs accretion rate

def load_accret(dir, width=None):
    import numpy as np
    
    acc = np.load(dir, allow_pickle=True).item()
    mass_cuts = acc['mass_cuts'][:-1]
    accret_med = acc['accret_med'] 
    z = np.round(acc['redshifts'], 3)

    if width is not None:
        if width =='std':
            acc_width = acc['accret_std']
        elif width == 'percentile':
            lowp_array  = acc['accret_low']
            highp_array = acc['accret_high']
            acc_width = highp_array - lowp_array
        return mass_cuts, z, accret_med, acc_width
    else:
        return mass_cuts, z, accret_med 

def plot_data(dir, snaps, acc, axs, cmap, norm, feat='depth'):
    import numpy as np
    import seaborn as sns
    from scipy.stats import pearsonr
    
    tot_x, tot_y = [], []
    
    acc_z, acc_mass_cuts, acc_data = acc[0], acc[1], acc[2]

    for snap in snaps:

        # Load data
        data = np.load(dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        z = np.round(data['z'], 3)
        print(f'snap {snap}: z = ', z)
        mass_cuts = data['mass_bins']
        print(mass_cuts)
        data = data['final_results']
        
        # In accretion rate, find the index corresponding to the current snap
        acc_snap_idx = np.where(acc_z == z)[0][0]
        print('acc z =', acc_z[acc_snap_idx], 'at idx', acc_snap_idx, acc_data.shape)
        
        # In accretion rate, find the index corresponding to the current mass cut
        comm_mass_cuts = np.intersect1d(acc_mass_cuts, mass_cuts)
        print(comm_mass_cuts)
        x_mass_idx = [np.where(acc_mass_cuts == m)[0][0] for m in comm_mass_cuts]
        print(x_mass_idx, acc_mass_cuts[x_mass_idx])
        y_mass_idx = [np.where(mass_cuts == m)[0][0] for m in comm_mass_cuts]
        print(y_mass_idx, mass_cuts[y_mass_idx])
        
        x = acc_data[acc_snap_idx, x_mass_idx]
        if feat == 'depth':
            y = data[y_mass_idx, 1, 1]
            y_min, y_max = data[y_mass_idx, 1, 0], data[y_mass_idx, 1, 2]
        elif feat == 'width':
            y = data[y_mass_idx, 2, 1]
            y_min, y_max = data[y_mass_idx, 2, 0], data[y_mass_idx, 2, 2]

        mask = x!=0
        x, y, y_min, y_max = x[mask], y[mask], y_min[mask], y_max[mask]
        
        # Sort according to x
        x = x[np.argsort(x)]
        y = y[np.argsort(x)]
        y_min = y_min[np.argsort(x)]
        y_max = y_max[np.argsort(x)]
        
        print('x:', x)
        print('y:', y)
        print('')
        
        axs.plot(x, y, color=cmap(norm(np.round(z, 3))), label=f'z = {z}')
        axs.errorbar(x, y, yerr=[y-y_min, y_max-y], color=cmap(norm(np.round(z, 3))), fmt='.')
        
        # Plot correlations
        # sns.set(style="whitegrid")
        # axs = sns.regplot(x=x, y=y, ci=95, scatter=False, line_kws={"color": cmap(norm(np.round(z, 3)))})
        
        tot_x = np.concatenate((tot_x, x))
        tot_y = np.concatenate((tot_y, y))
    
    return tot_x, tot_y

def plot_stats(dir, snaps, axs, cmap, norm, feat_idx):
    import numpy as np
    for snap in snaps:
        
        # Load data
        data = np.load(dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        z = data['z']
        
        mass_cuts = data['mass_bins']
        mass_cuts = [10**(10+m) for m in mass_cuts]
        
        data = data['final_results']

        axs.plot(mass_cuts, data[:, feat_idx, 1], 
                color=cmap(norm(np.round(z, 3))), alpha=0.2
                # label=f'z = {np.round(z, 1)}'
                )
        # axs.fill_between(mass_cuts, data[:, feat_idx, 0], data[:, feat_idx, 2], 
        #                  color=cmap(norm(np.round(z, 3))),
        #                  alpha=0.1)
        axs.errorbar(mass_cuts, data[:, feat_idx, 1],
                    yerr=[data[:, feat_idx, 1]-data[:, feat_idx, 0], 
                        data[:, feat_idx, 2]-data[:, feat_idx, 1]],
                    color=cmap(norm(np.round(z, 3))),
                    fmt='.')
    return axs