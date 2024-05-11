import numpy as np

radial_bins = np.logspace(0,np.log10(3*1),10)
radial_centers = 10**( (np.log10(radial_bins[1:])+np.log10(radial_bins[:-1]))/2 )

print(radial_bins)
print(radial_centers)