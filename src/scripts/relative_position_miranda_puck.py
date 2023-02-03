#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import paths
import pandas as pd
from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta

'''
Import the H-band detection image of Mab
Import the first frame, to which the offsets in all the other frames are tied
    Check which frame that is (they may or may not be time-sorted) using EXPSTART in the header
Put a dot on Mab's actual observed position in the frame
Pull the expected offsets of Mab, Puck, and Miranda, to make a triangle
See if that triangle matches the true position of Mab

Second test: flip the problem by overlaying the expected path of Miranda and Puck relative to 
    Mab on top of the detection image
'''

hdul = fits.open(paths.data / 'results/urh_Mab_2019-10-28_test3.fits')
single_frame_hdul = fits.open(paths.data / 'reduced/2019oct28/urh156.fits')
hdr = hdul[0].header
stack_data = hdul[0].data
frame_data = single_frame_hdul[0].data

obsdate = hdr["DATE-OBS"].strip(", \n")
obsdate = datetime.strftime(datetime.strptime(obsdate, "%Y-%m-%d"), "%Y-%b-%d")
start_time = obsdate + " " + hdr["EXPSTART"][:5]
obscode = '568'
date = '2019-10-28'
tstart = date+' 00:00'
tend = date+' 23:59'
pixscale_arcsec = 0.009971 #arcsec

# get positions of Puck, Miranda from the single frame and Mab from the stack
mab_position =  (193, 397) #(200, 234.5) #(274, 7.5)
miranda_position = (62, 481) #(49, 335) #(66, 109)
puck_position = (301, 738) #(242, 567) #(230, 281)
mab_motion = (-0.200581686891951, 6.719486510881566) #(-1.4040718082439412, 6.218032293651561) #(-3.1090161468258657, 4.8139604854076765) #from mab_path_on_detector.py

print(single_frame_hdul[0].header["EXPSTART"], start_time[-5:]) #these must agree


def extract_xy_from_horizons(id, start_times):
    '''
    id: Horizons object identifier
    start_times: must be list of times, each time has format "YYYY-Bbb-DD HH:MM"
    '''
    horizons_obj = Horizons(
        id=id,
        location=obscode,
        epochs={"start": tstart, "stop": tend, "step": "1m"},
    )
    ephem = horizons_obj.ephemerides(quantities=6).to_pandas()
    ephem = ephem.set_index(pd.DatetimeIndex(ephem["datetime_str"]))
    
    # match ephemeris time with time in fits header
    ephem_lines = ephem.loc[start_times]
    x_shift = ephem_lines["sat_X"].to_numpy() / pixscale_arcsec
    y_shift = ephem_lines["sat_Y"].to_numpy() / pixscale_arcsec
    
    return -x_shift, y_shift
    
# test extract multiple from Horizons
start_times = [start_time, start_time[:-1]+'9']
print(extract_xy_from_horizons('Puck', start_times))

## get ephemerides from Horizons. quantity 6 is the satellite relative position to parent in arcsec
x_shifts = []
y_shifts = []
for moon in ['Mab', 'Puck', 'Miranda']:
    
    x_shift, y_shift = extract_xy_from_horizons(moon, [start_time])
    x_shifts.append(x_shift)
    y_shifts.append(y_shift)
    
puck_mab_vector = [x_shifts[0] - x_shifts[1], y_shifts[0] - y_shifts[1]]
puck_miranda_vector = [x_shifts[2] - x_shifts[1], y_shifts[2] - y_shifts[1]]
miranda_mab_vector = [x_shifts[0] - x_shifts[2], y_shifts[0] - y_shifts[2]]

# plot the frame with these vectors on top
fig, ax = plt.subplots(1,1, figsize = (8,8))

ax.imshow(frame_data, origin='lower', vmin=0, vmax = 500, cmap='Greys_r')
ax.quiver(*puck_position, puck_miranda_vector[0], puck_miranda_vector[1], color = 'red',
            angles='xy', scale_units='xy', scale=1,
            label = 'Expected Relative Positions')
ax.quiver(*puck_position, puck_mab_vector[0], puck_mab_vector[1], color = 'red',
            angles='xy', scale_units='xy', scale=1,)
ax.quiver(*miranda_position, miranda_mab_vector[0], miranda_mab_vector[1], color = 'red',
            angles='xy', scale_units='xy', scale=1,)
            
ax.quiver(*mab_position, mab_motion[0], mab_motion[1], color = 'steelblue', label = 'Mab motion vector')

# true position of Mab
ax.scatter(mab_position[0], mab_position[1], color = 'steelblue', marker = 'o',
            label='Detected Position of Mab')

ax.set_ylim([350, 800])
ax.set_xlim([0, 400])
ax.set_xlabel('Pixels')
ax.set_ylabel('Pixels')
ax.legend()

fig.savefig(paths.figures / "expected_position_test3.png", dpi=300)
plt.show()
plt.close()



