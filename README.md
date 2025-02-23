# TNG

The codes in this branch are designed for simulation [IllustrisTNG](https://www.tng-project.org/) and [MillenniumTNG](https://www.mtng-project.org/).

## code

Install [illustris_python](https://github.com/illustristng/illustris_python) to ../code/illustris_python

### ../code/

* Adjust_profile_data_unit.py

  This script converts the data from TNG300 from comoving unit to phys unit, and the data from MTNG from Mpc to Kpc. This script should only be run once.

  - [ ] Double check the original MTNG DM halo profile data unit and decide where to store it. Check the code on mit branch.

### ../code/DMhalo_TNG300/

* 1_subset.py --boxsize --res --snapnum --mass_range --method --save_root_dir

  * one_halo_hist.py 

* check_halo_number.py --boxsize --res --mass_range

* 2_compile_all.py

  * compile_chunk_profile.py --snpnum --bin_start --bin_end

* 4_add_accretion.py

* 5_add_formation_time.py

### ../code/tree_TNG300/

The scripts in this folder should be run by the following order. All saved csv files have the same data structure and each entry in the csv files corresponds to the same halo. 

* table_mass.py --sim_type --bin_start --bin_end

  This script computes the DM halos masses history and saves as csv file.

  - [ ] How to find the info of mass bin start and bin end?

* table_snap_idx.py --sim_type
  
  This script computes the DM halos local index at snapshot X and the linking subhalo masses, and saves as csv files.

* table_accretion_new.py --sim_type

  This script computes the DM halos accretion rates history per dynamical time and saves as csv file.

### ../code/bootstrap/

* bootstrap_mass.py --sim --snapnum

* bootstrap_accret.py --sim --snapnum --accret_start --accret_end

  - [x] How to find the info of accretion rate bin start and bin end? see joob/TNG300/boots/

* boootstrap_z.py --sim --snapnum

### ../code/paper_plot/