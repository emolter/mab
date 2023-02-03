#!/usr/bin/env python

'''
TO DO:
- expand this to automatically choose proper observing blocks for each moon
    in order to maximize SNR of the detection
- expand this to automatically build Mab_full, Mab_block45, and Mab_block67
'''

from shift_stack_moons.shift_stack_moons import *
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from astroquery.jplhorizons import Horizons
import paths

# find filenames
data_dir = paths.data / 'reduced/2019oct28'
data_files = os.listdir(data_dir)
stem = 'urk'
obscode = '568'
code = 'Mab'
date = '2019-10-28'
tstart = date+' 00:00'
tend = date+' 23:59'
perturbation_experiment = False #set this flag to True to re-build mock data in paths.data / perturbation_experiment
pixscale = 0.009971

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
#select which observing blocks you want
#fnames = fnames_sorted[block_starts[5]+1:] #block_ends[-1]] #Mab is found in block_starts[5]+1 to end of block 7 for H-band
fnames = fnames_sorted

# for testing, make a different frame the zeroth one
#fnames = np.concatenate([fnames[12:], fnames[:12]])
#fnames = fnames[::-1]


## get ephemeris from Horizons. quantity 6 is the satellite relative position to parent in arcsec
horizons_obj = Horizons(
    id=code,
    location=obscode,
    epochs={"start": tstart, "stop": tend, "step": "1m"},
)
ephem = horizons_obj.ephemerides(quantities=6).to_pandas()
ephem = ephem.set_index(pd.DatetimeIndex(ephem["datetime_str"]))


if not perturbation_experiment:
    # do shift-and-stack and write
    fits_out = shift_and_stack(fnames, ephem, pixscale=0.009971, difference=True, edge_detect=True)
    fits_out.write(outfname)


if perturbation_experiment:
    # do a perturbation experiment to prove we don't see spurious sources
    for i in range(100):
        outfname2 = paths.data / f"perturbation_experiment/{code}_{date}_{i}.fits"
        fits_out = shift_and_stack(fnames, ephem, pixscale=pixscale, difference=True, edge_detect=False, perturbation_mode=True)
        fits_out.write(outfname2)
