#!/usr/bin/env python

'''
ready to ship
'''

from shift_stack_moons import shift_and_stack
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from astroquery.jplhorizons import Horizons
import paths


constants_dict = {
    'Mab':{
        'H':{'stem':'urh', 'reverse':True, 'canny':False, 'start_idx':41, 'end_idx':64}, #urh156
        'Kp':{'stem':'urk', 'reverse':True, 'canny':True, 'start_idx':0, 'end_idx':24}, #urk146
        },
    'Ophelia':{
        'H':{},
        'Kp':{},
        },
    'Cordelia':{
        'H':{},
        'Kp':{},
        },
    'Perdita':{
        'H':{},
        'Kp':{},
        },
}

# find filenames
data_dir = paths.data / 'reduced/2019oct28'
data_files = os.listdir(data_dir)
obscode = '568'
date = '2019-10-28'
tstart = date+' 00:00'
tend = date+' 23:59'
pixscale = 0.009971

np.random.seed(9)
if not os.path.exists(paths.data / "perturbation_experiment"):
    os.mkdir(paths.data / "perturbation_experiment")


for code in ['Mab']:
    ## get ephemeris from Horizons. quantity 6 is the satellite relative position to parent in arcsec
    horizons_obj = Horizons(
        id=code,
        location=obscode,
        epochs={"start": tstart, "stop": tend, "step": "1m"},
    )
    ephem = horizons_obj.ephemerides(quantities=6).to_pandas()
    ephem = ephem.set_index(pd.DatetimeIndex(ephem["datetime_str"]))
    
    for band in ['H', 'Kp']:
        print(f'Starting {band} band')
        stem = constants_dict[code][band]['stem']
        start_idx = constants_dict[code][band]['start_idx']
        end_idx = constants_dict[code][band]['end_idx']
        reverse = constants_dict[code][band]['reverse'] #do you want to input the files backward so the moon is stacked to the last position instead of the first position?
        canny = constants_dict[code][band]['canny']
        
        outfname = paths.data / f"results/{stem}_{code}_{date}.fits"
        
        fnames = [data_dir / s for s in data_files if s.startswith(stem)]
        
        # split into different observing blocks
        data_indices = np.array([int(s.split('.')[0][3:]) for s in data_files if s.startswith(stem)])
        sorti = np.argsort(data_indices)
        indices_sorted = np.sort(data_indices)
        diff = indices_sorted[1:] - indices_sorted[:-1]
        block_ends = np.argwhere(diff > 1).T[0]+1
        block_starts = np.concatenate([np.array([0]), block_ends])
        block_ends = np.concatenate([block_ends, np.array([len(fnames)])])
        fnames_sorted = np.array(fnames)[sorti]
        
        #select which observations you want - can separate them into blocks using block_starts and block_ends
        fnames = fnames_sorted[start_idx:end_idx] #block_ends[-1]] #Mab is found in block_starts[5]+1 to end of block 7 for H-band
        
        # for testing, make a different frame the zeroth one
        if reverse:
            fnames = fnames[::-1]
            
        for i in range(100):
            print(f"starting simulation {i} out of 100")
            outfname2 = paths.data / f"perturbation_experiment/{stem}_{code}_{date}_{i}.fits"
            fits_out = shift_and_stack(fnames, ephem, pixscale=pixscale, 
                        difference=True, edge_detect=canny, perturbation_mode=True, diagnostic_plots = False)
            fits_out.write(outfname2)