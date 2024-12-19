import numpy as np

def load_data(dir, snaps):
    
    for isnap, snap in enumerate(snaps): # from low z to high z (present)
            
        # Load data
        data = np.load(dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        z = data['z']
        
        mass_bins = data['mass_bins']
        mass_bins = [10**(10+m) for m in mass_bins]
        num_bins = len(mass_bins)
        
        data = data['final_results']

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

def load_z(dir, snaps):
    data = np.load(dir+f'/snap_{min(snaps)}_Rsp_stats.npy', allow_pickle=True).item()
    z_i = np.round(data['z'], 2)
    data = np.load(dir+f'/snap_{max(snaps)}_Rsp_stats.npy', allow_pickle=True).item()
    z_f = np.round(data['z'], 2)
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

def delta_c(z):
    import astropy.units as u
    from astropy.constants import M_sun, G
    from astropy.cosmology import FlatLambdaCDM
    
    # Define own cosmology
    h = 0.6774
    Om0 = 0.3089
    Ob0 = 0.0486
    cosmo = FlatLambdaCDM(H0 = h * 100 * u.km / u.s / u.Mpc, 
                            Om0=Om0, Ob0=Ob0, Tcmb0=2.725)
    
    # Compute the critical density
    H = cosmo.H(z) # km / (Mpc s)
    H = H.to(u.m / u.s / u.m)
    rho_c = 3 * H**2 / ( 8 * np.pi * G ) 
    rho_c = rho_c.to(M_sun / u.kpc**3) # Msun / kpc**3
    
    # Compute the characteristic overdensity
    delta_rho = rho_c * 1.686
    
    return delta_rho

def character_mass(char_r, char_rho):
    return 4/3 * np.pi * char_r**3 * char_rho