# TNG - TNG

The codes in this branch are designed for simulation [TNG](https://www.tng-project.org/).

## code

Install [illustris_python](https://github.com/illustristng/illustris_python) to ../code/illustris_python

### DMhalo

../code/DMhalo

* 1_subset.py --boxsize --res --snapnum --mass_range --method --save_root_dir

  * one_halo_hist.py 

* check_halo_number.py --boxsize --res --mass_range

* 2_compile_all.py

  * compile_chunk_profile.py --snpnum --bin_start --bin_end

### DMhalo_para

../code/DMhalo_para

* 1_subset.py --boxsize --res --snapnum --bin_start --bin_end --method --save_root_dir

  * one_halo_hist.py 

* 3_bootstrap_all.py

  * bootstrap.py --Nsample --Nboots --snapnum --bin_end

### stack_profiles (useless now)

* 2_stacks.py ---boxsize --res 

  * stacked_density_profiles.py --boxsize --res --snapnum --bin_start --bin_end --root_dir

  * profiles_time_evolution.py --boxsize --res --bin_start --bin_end 

  * Rsp_data_sort.py --boxsize --res --root_dir

### plot

* mass_hist.py

* Rspfeats_mass.py --feat_idx --bin_start --bin_end

* Rspfeats_redshift.py --feat_idx --bin_start --bin_end

* median_profiles.py