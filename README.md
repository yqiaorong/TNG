# TNG - MillenniumTNG

The codes in this branch are designed for simulation [MillenniumTNG](https://www.mtng-project.org/).

## code

Install [illustris_python](https://github.com/illustristng/illustris_python) to ../code/illustris_python

* add_accret.py

  This script adds the accretion rate to the corresponding DM halo density profile dataset.

### tree

../code/tree

The following scripts compute the accretion rates of DM halos.

* mass_table2.py 

  This script generates the mass table and the local index table of DM halos at each snapshot.

  e.g., Group_M_Mean200[index_table[i]] = mass_table[i]

* accretion_table2.py
  
  This script calculates the accretion rate of DM halos per dynamical time.