def compt_density_profile(coordinates, haloCM, halo_R_Mean200, dimension,
                          radius_range=[0.01, 5], number_of_bins=85):
    """This function computes the density profile of halos.
    
    INPUT:
    coordinates:    2D array with shape (N, 3)              [free unit]
    haloCM:         1D array with shape (3,)                [free unit]
    halo_R_Mean200: float
    dimension:      float                                   [ckpc^3/h]
    radius_range:   list with two fractional radius         [dimensionless]
    number_of_bins: float                       
    
    RETURN:
    density_bins:   1D array with shape (number of bins,)   [input unit^(-3)]
    radius_bins:    1D array with shape (number of bins,)   [input unit]
    """
    
    import numpy as np
    # radial_coordinates = np.linalg.norm(coordinates - haloCM, axis=1) 
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
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

# def compt_density_profile_hist(coordinates, mass_weights, haloPos, radial_bins):
    
#     radii = np.sqrt(np.sum(coordinates-haloPos, axis=1)**2) # [ckpc/h]
#     radial_volumes = 4/3*np.pi * (radial_bins[1:]**3 - radial_bins[:-1]**3) # [(ckpc/h)^3]
#     densities = np.histogram(radii,radial_bins,weights=mass_weights)[0] / radial_volumes # [Msun/h / (ckpc/h)^3]
#     return densities

import numpy as np
from scipy.optimize import curve_fit

def NSW_profile(r,
                rho_0, 
                R_s):
    return rho_0 / ((1 + (r/R_s)**2)*(r/R_s))

def fit_NSW_profile(bin_centers, densities, R_200_mean):
    
    def wrapped_NSW_profile(r:float, 
                        rho_0:float, 
                        R_s:float
                        ):
        return NSW_profile(r=r, rho_0=rho_0, R_s=R_s)
    
    
    base_p0 = (max(densities), R_200_mean)

    popt, pcov = curve_fit(NSW_profile,
                            bin_centers,
                            densities,
                            p0=base_p0,
                            maxfev=100000)
    perr = np.diag(pcov) ** 0.5
    
    def chi_square(p):
        return np.sum((p - densities)**2) / (len(p) - len(popt))

    # Did we actually get a good fit? If not, we should dump this bootstrapping.
    predicted_values = wrapped_NSW_profile(bin_centers, *popt)
    new_chi_square = chi_square(predicted_values)
    # If we get a worse fit after tuning, cancel this one
    old_chi_square = chi_square(wrapped_NSW_profile(bin_centers, *base_p0))
    # print(new_chi_square, old_chi_square)
    
    # if old_chi_square < new_chi_square:
    #     raise RuntimeError("Extremely poor fit for this bootstrap")
    # elif new_chi_square > 2.0:
    #     print("Bad Fit!")
    #     raise RuntimeError("Extremely poor fit for this bootstrap")

    return (
        popt,perr,
        )  