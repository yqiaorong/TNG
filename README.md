# TNG

The codes in this branch are designed for simulation [MillenniumTNG](https://www.mtng-project.org/).

# code

Install [illustris_python](https://github.com/illustristng/illustris_python) to ../code/illustris_python

### ../code/DMhalo_MTNG/

* 1_serial.py --sim --snapnum --start_chunk_idx --bin_start --bin_end

  * subset_tgt_chunk-gpu.py

* 2_compile_chunk_profile2.py

The MTNG density profiles data are saved in [data parent dir]/DMhalo_density_profiles/

They should be transferred to MPI branch and saved in TNG/DMhalo_density_profiles_phys/

Personal note: most of current data on mit branch is QUITE old so should NOT execute the above TRANSFER commands unless density profiles are ALL recomputed by rerunning above updated scripts! 