# TNG - MillenniumTNG

The codes in this branch are designed for simulation [MillenniumTNG](https://www.mtng-project.org/).

# code

Install [illustris_python](https://github.com/illustristng/illustris_python) to `../code/illustris_python`

### DMhalo_MTNG

`../code/DMhalo_MTNG/`

* 1_serial.py --sim --snapnum --start_chunk_idx --bin_start --bin_end

  * subset_tgt_chunk-gpu.py

* 2_compile_chunk_profile2.py --sim --snapnum --bin_start --bin_end

The table below gives the info of halo mass bins:

1. DM-Arepo

| snapnum | bin_start | bin_end |
|---------|-----------|---------|
| 129     | 3         | 4       |
| 151     | 3         | 4.5     |
| 179     | 3         | 4.5     |
| 214     | 3         | 5       |
| 237     | 3         | 5       |
| 264     | 3         | 5.5     |

2. Hydro-Arepo

| snapnum | bin_start | bin_end |
|---------|-----------|---------|
| 129     | 3         | 4       |
| 151     | 3         | 4.5     |
| 179     | 3         | 4.5     |
| 214     | 3         | 5       |
| 237     | 3         | 5       |
| 264     | 3         | 5.5     |

The MTNG density profiles data are saved in `[user name]/DMhalo_density_profiles/`

They should be transferred to MPI branch and saved in `/TNG/DMhalo_density_profiles_phys/`

> ⚠️ Personal note: most of current data on mit branch is QUITE old so should NOT execute the above TRANSFER commands unless density profiles are ALL recomputed by rerunning above updated scripts! 