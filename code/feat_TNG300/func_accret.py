def calc_tdyn(snap_dict, current_z):
    import numpy as np
    import astropy.units as u
    from astropy.cosmology import z_at_value, FlatLambdaCDM
    
    # Define own cosmology
    h = 0.6774
    Om0 = 0.3089
    Ob0 = 0.0486
    cosmo = FlatLambdaCDM(H0 = h * 100 * u.km / u.s / u.Mpc, 
                            Om0=Om0, Ob0=Ob0, Tcmb0=2.725)
    
    t = cosmo.age(current_z) # The cosmological time at current snapshot / z
    H = cosmo.H(current_z)   # The Hubble parameter at current snapshot / z
    Om = cosmo.Om(current_z) # The Omega m at current snapshot / z
    
    # Cosmological dynamic time
    t_H   = 1 / H
    t_dyn = t_H / (5 * np.sqrt(Om))
    t_dyn = t_dyn.to('yr').value * 1E-9 # t_dyn in unit Gyr
    print(f"One dynamical time before current z = {current_z} is {t_dyn} Gyr.")
    
    # Find one dynamical time before current z
    p_z = z_at_value(cosmo.age,  t - t_dyn * u.Gyr) # t_dyn in unit Gyr
    prev_snap = min(snap_dict, key=lambda snap: abs(snap_dict[snap] - p_z))
    prev_z = snap_dict[prev_snap]
    
    print('The previous redshift: ', prev_snap, prev_z)

    return prev_snap