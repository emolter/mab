#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib import cm
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

ready to ship!
'''

code = 'Mab'
band = 'H'

constants_dict = {
    'Mab':{
        'H':{'stem':'urh', 'frameid':156, 'xpos':187, 'ypos':395, 'Miranda_pos':np.array([62, 481]), 'Puck_pos':np.array([301, 738]), 'Mab_motion':np.array([-0.200581686891951, 6.719486510881566])}, #urh156
        'K':{'stem':'urk', 'frameid':146, 'xpos':188, 'ypos':307, 'Miranda_pos':np.array([54,403]), 'Puck_pos':np.array([267, 650]), 'Mab_motion':np.array([-0.9026175910138363, 6.518904823989544])}, #urk146
        },
    'Ophelia':{
        'H':{},
        'K':{},
        },
    'Cordelia':{
        'H':{},
        'K':{},
        },
    'Perdita':{
        'H':{},
        'K':{},
        },
}

obscode = '568'
date = '2019-10-28'
tstart = date+' 00:00'
tend = date+' 23:59'
pixscale_arcsec = 0.009971 #arcsec
stem = constants_dict[code][band]['stem']
frameid = constants_dict[code][band]['frameid']

hdul = fits.open(paths.data / f'results/{stem}_{code}_2019-10-28.fits')
single_frame_hdul = fits.open(paths.data / f'reduced/2019oct28/{stem}{frameid}.fits') #start is urh104, end is urh156
hdr = hdul[0].header
stack_data = hdul[0].data
frame_data = single_frame_hdul[0].data
obsdate = hdr["DATE-OBS"].strip(", \n")
obsdate = datetime.strftime(datetime.strptime(obsdate, "%Y-%m-%d"), "%Y-%b-%d")
start_time = obsdate + " " + hdr["EXPSTART"][:5]

# get positions of Puck, Miranda from the single frame and Mab from the stack
mab_position = np.array([constants_dict[code][band]['xpos'], 
                        constants_dict[code][band]['ypos']])
miranda_position = constants_dict[code][band]['Miranda_pos']
puck_position = constants_dict[code][band]['Puck_pos'] 
mab_motion = constants_dict[code][band]['Mab_motion']  #(-1.4040718082439412, 6.218032293651561) #(-3.1090161468258657, 4.8139604854076765) #from mab_path_on_detector.py

mab_start_position = mab_position #- mab_motion/2

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
    

## get ephemerides from Horizons. quantity 6 is the satellite relative position to parent in arcsec
x_shifts = []
y_shifts = []
for moon in [code, 'Puck', 'Miranda']:
    
    x_shift, y_shift = extract_xy_from_horizons(moon, [start_time])
    x_shifts.append(x_shift)
    y_shifts.append(y_shift)
    
puck_mab_vector = np.array([x_shifts[0] - x_shifts[1], y_shifts[0] - y_shifts[1]])[:,0]
puck_miranda_vector = np.array([x_shifts[2] - x_shifts[1], y_shifts[2] - y_shifts[1]])[:,0]
miranda_mab_vector =  np.array([x_shifts[0] - x_shifts[2], y_shifts[0] - y_shifts[2]])[:,0]

expected_mab_position_1 = puck_position + puck_mab_vector
expected_mab_position_2 = miranda_position + miranda_mab_vector
#print(expected_mab_position_1, expected_mab_position_2, mab_position)

# plot the frame with these vectors on top
fig, ax = plt.subplots(1,1, figsize = (8,8))

offset = 100
data_show = frame_data+offset
data_show[data_show < offset] = offset
ax.imshow(data_show, origin='lower', norm=LogNorm(vmin=-1 + offset, vmax=1e3 + offset, clip=False), cmap='Greys_r')
ax.quiver(*puck_position, puck_miranda_vector[0], puck_miranda_vector[1], color = 'red',
            angles='xy', scale_units='xy', scale=1,
            label = 'Expected Relative Positions')
ax.quiver(*puck_position, puck_mab_vector[0], puck_mab_vector[1], color = 'red',
            angles='xy', scale_units='xy', scale=1,)
ax.quiver(*miranda_position, miranda_mab_vector[0], miranda_mab_vector[1], color = 'red',
            angles='xy', scale_units='xy', scale=1,)
            
ax.quiver(*mab_position, mab_motion[0], mab_motion[1], color = 'steelblue', label = f'{code} motion vector')

## true position of Mab
ax.scatter(mab_start_position[0], mab_start_position[1], 
            edgecolor = 'magenta', marker = 'o', facecolor='none', s=50,
            label=f'Detected Position of {code}')

ax.set_ylim([350, 800])
ax.set_xlim([0, 400])
ax.set_xlabel('Pixels')
ax.set_ylabel('Pixels')
ax.legend()

fig.savefig(paths.figures / f"expected_position_{stem}_{code}_oneframe.png", dpi=300)
#plt.show()
plt.close()

## do the flipped problem here: assume Mab position is correct and static on the detector
##   then where should the paths of Miranda and Puck go?
time0 = '2019-10-28 10:36'
time1 = '2019-10-28 12:48'
alltimes = pd.date_range(time0, time1, 133)

miranda_eph = extract_xy_from_horizons('Miranda', alltimes)
puck_eph = extract_xy_from_horizons('Puck', alltimes)
mab_eph = extract_xy_from_horizons(code, alltimes)

puck_vec = (puck_eph[0] - mab_eph[0], puck_eph[1] - mab_eph[1])
miranda_vec = (miranda_eph[0] - mab_eph[0], miranda_eph[1] - mab_eph[1])

origin = np.ones((miranda_vec[0].shape[0],2)) * np.array([mab_position[0], mab_position[1]])
fig, ax = plt.subplots(1,1, figsize = (8,8))

cmap = cm.get_cmap("Greys_r").copy()
cmap.set_under('k')
cmap.set_over('cyan')
offset = 100
data_show = stack_data+offset
data_show[data_show < offset] = offset
cim = ax.imshow(data_show, origin='lower', cmap=cmap, norm=LogNorm(vmin=-1 + offset, vmax=7e3 + offset, clip=False))

#ax_divider = make_axes_locatable(ax)
#cax = ax_divider.append_axes("right", size="7%", pad="2%")
#tick_labels = np.array([0, 300, 1000, 5000])
#cbar = fig.colorbar(cim, cax=cax, label=r'Total Flux (cts)', ticks = tick_labels + offset)
#cbar.ax.set_yticklabels(tick_labels)

#ax.quiver(*origin.T, puck_vec[0], puck_vec[1], color = 'red',
#            angles='xy', scale_units='xy', scale=1,
#            label = 'Expected Relative Positions')
#ax.quiver(*origin.T, miranda_vec[0], miranda_vec[1], color = 'red',
#            angles='xy', scale_units='xy', scale=1,)
ax.scatter(puck_vec[0] + origin[:,0], puck_vec[1] + origin[:,1], 
            color = 'red', marker = '.', s = 10,
            label='Expected Trajectory of Puck')
ax.scatter(miranda_vec[0] + origin[:,0], miranda_vec[1] + origin[:,1], 
            color = 'orange', marker = '.', s=10,
            label='Expected Trajectory of Miranda')

# true position of Mab
ax.scatter(mab_start_position[0], mab_start_position[1], 
            edgecolor = 'magenta', facecolor='none', marker = 'o', s = 50,
            label=f'Detected Position of {code}')

# hack to make Miranda brightest spots show up on 
ax.scatter(0,0, marker='o', color='cyan', label='Brightest Pixels in Miranda')

ax.set_ylim([350, 800])
ax.set_xlim([0, 400])
ax.set_xlabel('Pixels')
ax.set_ylabel('Pixels')
ax.legend()

fig.savefig(paths.figures / f"expected_position_{stem}_{code}.png", dpi=300)
#plt.show()
plt.close()




