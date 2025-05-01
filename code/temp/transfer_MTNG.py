import os
import numpy as np

snaps = [#129, 151, 179, 214, 237, 
         264]

for snap in snaps:
    from_dir = f'buffer/MTNG-L500-4320-A/snap_{snap}/final_densities/'
    from_fname = os.listdir(from_dir)[0]
    print(f'From fname: {from_fname}')
    from_data = np.load(from_dir + from_fname, allow_pickle=True).item()
    print(from_data.keys())
    
    to_dir = f'result/DMhalo_density_profiles/MTNG/Hydro-Arepo/MTNG-L500-4320-A/snap_{snap}/final_densities/'
    to_fname = os.listdir(to_dir)[0]
    print(f'To fname: {to_fname}')
    to_data = np.load(to_dir + to_fname, allow_pickle=True).item()
    print(to_data.keys())
    
    
    # # Transfer accretions
    # to_data['accretions'] = from_data['accretions']
    # # Transfer GroupNum
    # to_data['GroupNum'] = from_data['GroupNum']
    # # Transfer FirstSub
    # to_data['FirstSub'] = from_data['FirstSub']
    # # Transfer formzOLD
    # to_data['formzOLD'] = from_data['formzOLD']
    # # Transfer formzSub
    # to_data['formzSub'] = from_data['formzSub']
    # Transfer mergerz
    to_data['mergerz'] = from_data['mergerz']
    
    print(to_data.keys())
    
    # Save the modifed to data
    np.save(to_dir + to_fname, to_data)
    print('Saved!')
    print('')
    