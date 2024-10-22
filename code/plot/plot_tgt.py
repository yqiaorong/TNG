import os
import numpy as np

for idx in range(3):
    os.system(f'python3 code/plot/plot_tgt_mass.py --DM DM --feat_idx {idx}')
    os.system(f'python3 code/plot/plot_tgt_mass.py --DM Hydro --feat_idx {idx}')
    os.system(f'python3 code/plot/plot_tgt_redshift.py --DM DM --feat_idx {idx}')
    os.system(f'python3 code/plot/plot_tgt_redshift.py --DM Hydro --feat_idx {idx}')
    

cuts = np.linspace(1, 4, 7)
for cut in cuts:
    os.system(f'python3 code/plot/median_profiles_evolve.py --sim TNG300 --DM DM --mass_cut {cut}')
    
cuts = np.linspace(1.5, 4, 6)
for cut in cuts:
    os.system(f'python3 code/plot/median_profiles_evolve.py --sim TNG300 --DM Hydro --mass_cut {cut}')

cuts = np.linspace(3, 5, 6)
for cut in cuts:
    os.system(f'python3 code/plot/median_profiles_evolve.py --sim MTNG --DM DM --mass_cut {cut}')
    
cuts = np.linspace(3, 5, 5)
for cut in cuts:
    os.system(f'python3 code/plot/median_profiles_evolve.py --sim MTNG --DM Hydro --mass_cut {cut}')
    
    
os.system('python3 code/plot/plot_depth_vs_accret.py --DM Hydro')
os.system('python3 code/plot/plot_depth_vs_accret.py --DM DM')


os.system('python3 code/plot/plot_width_vs_accret_width.py --DM DM --width percentile')
os.system('python3 code/plot/plot_width_vs_accret_width.py --DM DM --width std')
os.system('python3 code/plot/plot_width_vs_accret_width.py --DM Hydro --width percentile')
os.system('python3 code/plot/plot_width_vs_accret_width.py --DM Hydro --width std')