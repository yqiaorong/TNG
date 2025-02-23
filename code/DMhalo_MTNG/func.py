def compt_density_profile(coordinates, haloCM, halo_R_Mean200, dimension,
                          radius_range=[0.01, 5], number_of_bins=85):
    """This function computes the density profile of halos.
    
    INPUT:
    coordinates:    2D array with shape (N, 3)              [free unit]
    haloCM:         1D array with shape (3,)                [free unit]
    halo_R_Mean200: float                                   [free unit]
    dimension:      float                                   [ckpc^3/h]
    radius_range:   list with two fractional radius         [dimensionless]
    number_of_bins: float                       
    
    RETURN:
    density_bins:   1D array with shape (number of bins,)   [input unit^(-3)]
    radius_bins:    1D array with shape (number of bins,)   [input unit]
    """
    
    import numpy as np 
    radial_coordinates = distance(haloCM, coordinates, dimension)
    
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

def distance(x0, x1, dimensions):
    '''
    find distance between points in a periodic box
    x0,x1=coordinates or set of coordinates
    delta = distance between two points assuming no periodicity
    dimensions = box size *7500 ckpc/h*
    returns delta if delta<box size, else returns 1-delta as the distance between points
    '''
    import numpy as np
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))