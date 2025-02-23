# TNG - MillenniumTNG

The codes in this branch are designed for simulation [MillenniumTNG](https://www.mtng-project.org/).

## code

Install [illustris_python](https://github.com/illustristng/illustris_python) to `../code/illustris_python`

`../code/`

* add_accret.py --sim_type --snapnum

  This script adds the accretion rate to the corresponding DM halo density profile dataset.

* add_formation_time.py --sim --snapnum

  This script adds the formation time to the corresponding DM halo density profile dataset.

### tree_MTNG

`../code/tree_MTNG/`

* mass_table2.py --sim --restart

  This script generates the mass table and the local index table of DM halos at each snapshot.

  e.g., Group_M_Mean200[index_table[i]] = mass_table[i]

* accretion_table2.py --sim
  
  This script calculates the accretion rate of DM halos per dynamical time.