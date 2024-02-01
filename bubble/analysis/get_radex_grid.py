from spectralradex import radex
import numpy as np
from astropy import units as u 
from astropy import constants as c
from astropy.table import Table
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

import warnings
warnings.filterwarnings('ignore')

params = radex.get_example_grid_parameters()

params['tkin'] = np.arange(30, 150, 5)
params['h2'] = 10**np.arange(1, 8, 0.05)
params['linewidth'] = [1, 5, 10, 20]

N_h2 = 8.5e22 
X_hnco = 10**np.arange(-10, -5, 0.05)
params['cdmol'] = X_hnco * N_h2

params['molfile'] = 'hnco.dat'
params['fmin'] = 87.925
params['fmax'] = 87.926

renrun = True
if renrun:
    pool=Pool(36)
    grid_df = radex.run_grid(params, target_value="T_R (K)", pool=pool)
    grid_df.to_csv('./../data/radex/grid_df.csv')
else: 
    print('[INFO] No rerun')

grid = Table.read('./../data/radex/grid_df.csv')

grid['col0'].name = 'id'
grid['h2'].name = 'n_h2'
grid['cdmol'].name = 'N_hnco'
grid['(4_0_4)-(3_0_3)[87.92523962 GHz]'].name = 'I_hnco'

# Some the values did not converge are are missing, which makes assigning an Xhnco a little difficult
X_hnco_column = []
for id in tqdm(range(len(grid))):
    where = np.float32(grid['N_hnco'][id]) == np.float32(np.round(N_hnco))
    X_hnco_column += [np.round(X_hnco[np.where(where)[0][0]], 15)]
grid.add_column(X_hnco_column, name='X_hnco')

mask = (grid['I_hnco']>1)&(grid['I_hnco']<5)
grid.add_column(mask*1, name='mask')