"""This script uses the subhalo first progenitor information in group catalog to track the mass 
evolution of the first progenitor of the halos."""

import os
import h5py
import argparse
import numpy as np
import illustris_python as il
from tqdm import tqdm

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default=None,  type=str)
parser.add_argument('--restart', default=False, type=bool)
args = parser.parse_args()

print('')
print(f'>>> Mass table 2 <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



##############################################################################################
# I customise some functions
##############################################################################################

import six
from os.path import isfile,expanduser

def FProgPath(basePath, snapNum, chunkNum=0):
    """ Return absolute path to a group catalog HDF5 file (modify as needed). """
    gcPath = basePath + '/groups_%03d/' % snapNum
    filePath1 = gcPath + 'groups_%03d.%d.hdf5' % (snapNum, chunkNum)
    filePath2 = gcPath + 'subhalo_prog_%03d.%d.hdf5' % (snapNum, chunkNum)

    if isfile(expanduser(filePath1)):
        return filePath1
    return filePath2

def load_FProg_Objects(basePath, snapNum, gName, nName, fields):
    """ Load either halo or subhalo information from the group catalog. """
    result = {}

    # make sure fields is not a single element
    if isinstance(fields, six.string_types):
        fields = [fields]

    # load header from first chunk
    with h5py.File(FProgPath(basePath, snapNum), 'r') as f:

        header = dict(f['Header'].attrs.items())

        if 'N'+nName+'_Total' not in header and nName == 'subgroups':
            nName = 'subhalos' # alternate convention

        result['count'] = f['Header'].attrs['N' + nName + '_Total']

        if not result['count']:
            print('warning: zero groups, empty return (snap=' + str(snapNum) + ').')
            return result

        # if fields not specified, load everything
        if not fields:
            fields = list(f[gName].keys())

        for field in fields:
            # verify existence
            if field not in f[gName].keys():
                raise Exception("Group catalog does not have requested field [" + field + "]!")

            # replace local length with global
            shape = list(f[gName][field].shape)
            shape[0] = result['count']

            # allocate within return dict
            result[field] = np.zeros(shape, dtype=f[gName][field].dtype)

    # loop over chunks
    wOffset = 0

    for i in range(header['NumFiles']):
        f = h5py.File(FProgPath(basePath, snapNum, i), 'r')

        if not f['Header'].attrs['N'+nName+'_ThisFile']:
            continue  # empty file chunk

        # loop over each requested field
        for field in fields:
            if field not in f[gName].keys():
                raise Exception("Group catalog does not have requested field [" + field + "]!")

            # shape and type
            shape = f[gName][field].shape

            # read data local to the current file
            if len(shape) == 1:
                result[field][wOffset:wOffset+shape[0]] = f[gName][field][0:shape[0]]
            else:
                result[field][wOffset:wOffset+shape[0], :] = f[gName][field][0:shape[0], :]

        wOffset += shape[0]
        f.close()

    # only a single field? then return the array instead of a single item dict
    if len(fields) == 1:
        return result[fields[0]]

    return result

def load_FProg_Subhalos(basePath, snapNum, fields=None):
    """ Load all subhalo information from the entire group catalog for one snapshot
       (optionally restrict to a subset given by fields). """

    return load_FProg_Objects(basePath, snapNum, "Subhalo", "subgroups", fields)

##############################################################################################
# Customisation ends here
##############################################################################################



basePath = f'/virgotng/mpa/MTNG/{args.sim}/output/'

# Save directory
save_dir = f'result/DMhalo_mass_table_new/{args.sim}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

Subhalo_fields = ['SubhaloGroupNr', 'SubhaloMass'] # A range from <0> to <total subhalos number> per current snapshot, RESET per snapshot!!!
Halo_fields = ['Group_M_Mean200', 'GroupFirstSub']



# Get the snapshot number of the last saved file
save_list = os.listdir(save_dir)
save_list = [item for item in save_list if item.endswith('_FPGr.npy')]
if len(save_list) == 0 or args.restart == True:
    min_snap = 264
    print(f'load halos in snap {min_snap}...')
    GroupFirstSub = il.groupcat.loadHalos(basePath, min_snap, fields='GroupFirstSub')
    SubhaloMass = il.groupcat.loadSubhalos(basePath, min_snap, fields='SubhaloMass')
    np.save(f'{save_dir}/snap_{min_snap}_FPGr.npy', range(GroupFirstSub.shape[0]))
    np.save(f'{save_dir}/snap_{min_snap}_SubMass.npy', SubhaloMass)
    del GroupFirstSub
    del SubhaloMass
else:
    snaps = [int(item.split('_')[1]) for item in save_list]
    min_snap = min([item for item in snaps])
    print(snaps, min_snap)
    # Remove the files that are not the latest snapshot
    for snap in snaps:
        if snap != min_snap and snap not in [51, 69, 94, 129, 151, 179, 214, 237, 264]:
            print(f'remove {snap}')
            os.remove(f'{save_dir}/snap_{snap}_FPGr.npy')
            os.remove(f'{save_dir}/snap_{snap}_SubMass.npy')
    print('')



# Start from the latest snapshot
print(f'Last saved snapshot: {min_snap}')
for snap in tqdm(range(min_snap, 15, -1)):
        
    ### Find the group index in the previous snapshot ###
    # -------------------------------------------------------------------------------------------------------
    if snap == min_snap:
        FPGr_indices = np.load(f'{save_dir}/snap_{min_snap}_FPGr.npy') 
        link_SubMass = np.load(f'{save_dir}/snap_{min_snap}_SubMass.npy')
    else:
        # Load the subhalo catalog in current snapshot
        print(f'load subhalos in snap {snap}...')
        Subhalos = il.groupcat.loadSubhalos(basePath, snap, fields=Subhalo_fields)
        SubhaloGroupNr = Subhalos['SubhaloGroupNr'] # Index into halos 
        SubhaloMass    = Subhalos['SubhaloMass']
        del Subhalos
        
        # Create the mapping
        SubGr_map = {Sub_idx: GrNr for Sub_idx, GrNr in enumerate(SubhaloGroupNr)} 
        SubMass_map = {Sub_idx: mass for Sub_idx, mass in enumerate(SubhaloMass)}
        
        # Map
        FPGr_indices = np.array([SubGr_map.get(Sub_idx, -1) for Sub_idx in FPSub_indices]) # Index into selected halos 
        link_SubMass = np.array([SubMass_map.get(Sub_idx, -1) for Sub_idx in FPSub_indices])
        np.save(f'{save_dir}/snap_{snap}_FPGr.npy', FPGr_indices) # This saved the group indices in the snapshot in filename
        np.save(f'{save_dir}/snap_{snap}_SubMass.npy', link_SubMass)
        
        # Remove the useless files
        if os.path.exists(f'{save_dir}/snap_{snap+1}_FPGr.npy') and snap+1 not in [51, 69, 94, 129, 151, 179, 214, 237, 264]:
            os.remove(f'{save_dir}/snap_{snap+1}_FPGr.npy')
        if os.path.exists(f'{save_dir}/snap_{snap+1}_SubMass.npy') and snap+1 not in [51, 69, 94, 129, 151, 179, 214, 237, 264]:
            os.remove(f'{save_dir}/snap_{snap+1}_SubMass.npy')
            
    print(FPGr_indices.shape)
    print(link_SubMass.shape)
    print('')
        
    ### Load the group catalog in current snapshot ###
    # -------------------------------------------------------------------------------------------------------
    print(f'load halos in snap {snap}...')
    Halos = il.groupcat.loadHalos(basePath, snap, fields=Halo_fields)
    GroupFirstSub   = Halos['GroupFirstSub']  
    Group_M_Mean200 = Halos['Group_M_Mean200']
    del Halos
    
    # Create the mapping
    GrMass_map = {GrNr: mass for GrNr, mass in enumerate(Group_M_Mean200)}
    
    if snap != min_snap:
        print(f'match? ({min(SubhaloGroupNr)}, {max(SubhaloGroupNr)}) < {GroupFirstSub.shape}')
        del SubhaloGroupNr
    
    # Map
    FPGrMass = np.array([GrMass_map.get(GrNr, -1) for GrNr in FPGr_indices])
    if snap in [51, 69, 94, 129, 151, 179, 214, 237, 264]:
       np.save(f'{save_dir}/snap_{snap}_FPGrMass.npy', FPGrMass)
       print(f'snap_{snap}_FPGrMass.npy saved!')
       
    
    ### Load the subhalo index in the previous snapshot by loading the subhalo catalog in current snaphot ###
    # -------------------------------------------------------------------------------------------------------
    print(f'load FP subhalos saved in snap {snap}...')
    FirstProgSubhaloNr = load_FProg_Subhalos(basePath, snap, fields='FirstProgSubhaloNr') # Index into subhalos in previous snapshot
    
    # Create the mapping
    FPSub_map = {Sub_idx: FPSub_idx for Sub_idx, FPSub_idx in enumerate(FirstProgSubhaloNr)} 
    del FirstProgSubhaloNr
    
    # Map
    FPSub_indices = np.array([FPSub_map.get(Sub_idx, -1) for Sub_idx in GroupFirstSub[FPGr_indices]]) # Index into selected subhalos in previous snapshot
    print(FPSub_indices.shape)