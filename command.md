python3 code/plot/mass_hist.py --sim DM-Arepo/MTNG-L500-4320-A/output --snap 129

mpirun -np 20 python3 code/DMhalo_para/1_subset_tgt_all.py --sim Hydro-Arepo/MTNG-L500-4320-A/output --bin_start 3.5 --bin_end 4 --snapnum 214

python3 code/DMhalo_tgt/1_subset_tgt_all.py --sim DM-Arepo/MTNG-L500-4320-A/output --bin_start 3 --bin_end 3.5 --snapnum 264

python3 code/DMhalo_tgt/2_compile_chunk_profile.py --sim Hydro-Arepo/MTNG-L500-4320-A/output --bin_start 3.5 --bin_end 4 --snapnum 214

python3 code/DMhalo_tgt/bootstrap_tgt.py --sim Hydro-Arepo/MTNG-L500-4320-A/output --bin_start 3.5 --bin_end 4 --snapnum 129

python3 code/plot/Rspfeats_mass.py --sim DM-Arepo/MTNG-L500-4320-A/output
python3 code/plot/Rspfeats_redshift.py --sim DM-Arepo/MTNG-L500-4320-A/output

/virgotng/mpa/MTNG/

DM-Arepo/MTNG-L500-4320-A/output
Hydro-Arepo/MTNG-L500-4320-A/output