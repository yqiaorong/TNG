# TNG

The codes in this branch are designed for simulation [IllustrisTNG](https://www.tng-project.org/) and [MillenniumTNG](https://www.mtng-project.org/).

## code

Install [illustris_python](https://github.com/illustristng/illustris_python) to `../code/illustris_python`

### DMhalo_TNG300

`../code/DMhalo_TNG300/`

* 1_subset.py --sim_type --snapnum --bin_start --bin_end

  * one_halo_hist.py 

* check_halo_number.py --sim_type --bin_start --bin_end

* 2_compile_all.py

  * compile_chunk_profile.py

- [ ] The above scripts have problems of finding the info of bin_start and bin_end as well.

* 3_comoving_to_phys.py

  This script converts the data from TNG300 from comoving unit to phys unit, and it should only be run once.

* 4_add_accretion.py --sim_type --snapnum

* 5_add_formation_time.py --sim_type --snapnum

### tree_TNG300

`../code/tree_TNG300/`

The scripts in this folder should be run by the following order. All saved csv files have the same data structure and each entry in the csv files corresponds to the same halo. 

* table_mass.py --sim_type --bin_start --bin_end

  This script computes the DM halos masses history and saves as csv file.

  - [ ] How to find the info of mass bin start and bin end?

* table_snap_idx.py --sim_type
  
  This script computes the DM halos local index at snapshot X and the linking subhalo masses, and saves as csv files.

* table_accretion_new.py --sim_type

  This script computes the DM halos accretion rates history per dynamical time and saves as csv file.

### bootstrap

`../code/bootstrap/`

* bootstrap_mass.py --sim --snapnum

* bootstrap_accret.py --sim --snapnum --accret_start --accret_end

  - [x] How to find the info of accretion rate bin start and bin end? see the tables below!

* boootstrap_z.py --sim --snapnum

The tables below summarise the range of accretion rates used in boostrapping:

1. TNG300/sim_205_1250_DM

  | snapnum | accret_start | accret_end |
  |---------|--------------|------------|
  | 8       | 0            | 0          |
  | 13      | 1            | 13         |
  | 17      | 1            | 13         |
  | 21      | 1            | 12         |
  | 25      | 0            | 14         |
  | 33      | 0            | 13         |
  | 40      | 0            | 14         |
  | 50      | 0            | 14         |
  | 67      | 0            | 13         |
  | 78      | 0            | 13         |
  | 99      | 0            | 7          |

2. TNG300/sim_205_1250_Hydro

  | snapnum | accret_start | accret_end |
  |---------|--------------|------------|
  | 8       | 0            | 0          |
  | 13      | 1            | 13         |
  | 17      | 1            | 13         |
  | 21      | 1            | 12         |
  | 25      | 1            | 13         |
  | 33      | 0            | 14         |
  | 40      | 0            | 14         |
  | 50      | 0            | 14         |
  | 67      | 0            | 13         |
  | 78      | 0            | 13         |
  | 99      | 0            | 7          |

3. MTNG/DM-Arepo

4. MTNG/Hydro-Arepo

### paper_plot

`../code/paper_plot/`