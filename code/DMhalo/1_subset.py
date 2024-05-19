import illustris_python as il
import os
import argparse
import numpy as np

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,   type=int)
parser.add_argument('--res',default=1250,      type=int)
parser.add_argument('--snapnum',default=99,    type=int)
parser.add_argument('--mass_range',default=3.5,type=float) # [10^{10+x} Msun/h]

parser.add_argument('--method', default=None,  type=str)
parser.add_argument('--save_root_dir',default='DMhalo_density_profiles',type=str)
args = parser.parse_args()

print('')
print(f'>>> DM halo density profiles subset <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Specify the snapshot
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(args.boxsize,args.res)
snapnum = args.snapnum
boxsize = args.boxsize
res     = args.res



# Save root dir
save_root_dir = args.save_root_dir +'_'+args.method
    
    

# Select a subset of DM halos from groupcat
group_fields = ['GroupPos', 'Group_M_Mean200', 'Group_R_Mean200']
Halos = il.groupcat.loadHalos(basePath, snapnum, fields=group_fields)

GroupPos        = Halos['GroupPos']        # [ckpc/h]
Group_M_Mean200 = Halos['Group_M_Mean200'] # [10^10 MSun/h]
Group_R_Mean200 = Halos['Group_R_Mean200'] # [ckpc/h]

# Using physical mass to select subset
subset_idx = np.where((Group_M_Mean200 >= 10**args.mass_range) & 
                      (Group_M_Mean200 < 10**(args.mass_range+0.5)))[0]
Ngroups_subset = subset_idx.shape[0]
print(f'In total, {Ngroups_subset} DM halos with mass 10^{args.mass_range+10} ~ 10^{args.mass_range+10.5} MSun in at snap {snapnum}')



# Iterate over DM halos
for i, idx in enumerate(subset_idx):
    if not os.path.exists(f'result/{save_root_dir}/sim_{boxsize}_{res}/snap_{snapnum}/densities/halo_{idx}.npy'):
        # Round values 
        x, y, z = np.round(GroupPos[idx, 0].item(), 0), np.round(GroupPos[idx, 1].item(), 0), np.round(GroupPos[idx, 2].item(), 0)
        R = np.round(Group_R_Mean200[idx].item(), 0)
        # Run the script
        os.system(f'python3 code/DMhalo/one_halo_hist.py'+
            f' --boxsize {boxsize} --res {res} --snapnum {snapnum} --groupnum {idx}'+
            f' --x {x} --y {y} --z {z}'+
            f' --M {Group_M_Mean200[idx]} --R {R}'+
            f' --save_root_dir {save_root_dir} --method {args.method}')
    else:
        print(f'At snap {snapnum}, DM halo local index {i+1}/{Ngroups_subset} already exists.')

print(f'All DM halos in the subset at snap {snapnum} are finished.')