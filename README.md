# TNG - MillenniumTNG

The codes in this branch are designed for simulation [MillenniumTNG](https://www.mtng-project.org/).

## code

Install [illustris_python](https://github.com/illustristng/illustris_python) to ../code/illustris_python

### DMhalo

../code/DMhalo

The following scripts compute DM halo densities in old (slow) way. 

* 1_subset.py --sim --snapnum --mass_range --save_root_dir

  * one_halo_hist.py 

* 3_bootstrap_all.py

  * bootstrap.py --sim --snapnum --Nsample

### DMhalo_tgt

The following scripts compute DM halo densities in new (fast) way. 

* 1_subset_tgt_all.py

  * subset_tgt_chunk.py

* 2_compile_chunk_profile.py

* bootstrap_tgt.py

### plot

../code/plot

* mass_hist.py --sim --snapnum --bin_start --bin_end 

* Rspfeats_mass.py

* Rspfeats_redshift.py