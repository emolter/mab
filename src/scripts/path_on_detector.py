#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import paths
import os
import pandas as pd
from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta
from pylanetary.navigation import *
from pylanetary.utils import Body
import astropy.units as u

'''
imshow a single detector frame, then overplot the path Mab took
    over the course of the observations
    one point per actual exposure
determine the average Mab vector over a single exposure using data from Horizons
compare that with the elongation of Mab

ready to ship!
'''

# find filenames
data_dir = paths.data / 'reduced/2019oct28'
data_files = os.listdir(data_dir)
obscode = '568'
code = 'Mab'
date = '2019-10-28'
tstart = date+' 00:00'
tend = date+' 23:59'
pixscale_arcsec = 0.009971 #arcsec
req = 25559 #km
rpol = 24973 #km

constants_dict = {
    'Mab':{'H':{'stem':'urh', 'canny':False, 'start_idx':41, 'end_idx':64}, #urh156
            'K':{'stem':'urk', 'canny':True, 'start_idx':0, 'end_idx':24}}, #urk146
    'Perdita':{'H':{'stem':'urh', 'canny':False, 'start_idx':0, 'end_idx':48}, 
            'K':{'stem':'urk', 'canny':True, 'start_idx':0, 'end_idx':11}}, 
    'Cupid':{'H':{'stem':'urh', 'canny':False, 'start_idx':0, 'end_idx':64}, 
            'K':{'stem':'urk', 'canny':True, 'start_idx':0, 'end_idx':24}}, 
        }
# Cupid on detector in all K-band frames but should be reversed
# Cupid on detector for frames 0-29, 48-64 in H-band and also should be reversed

## get ephemeris from Horizons. quantity 6 is the satellite relative position to parent in arcsec
horizons_obj = Horizons(
    id=code,
    location=obscode,
    epochs={"start": tstart, "stop": tend, "step": "1m"},
)
ephem = horizons_obj.ephemerides(quantities=6).to_pandas()
ephem = ephem.set_index(pd.DatetimeIndex(ephem["datetime_str"]))

exposure_starts = []
filters = ['K', 'H']
for filt in filters:
    stem = constants_dict[code][filt]['stem']
    start_idx = constants_dict[code][filt]['start_idx']
    end_idx = constants_dict[code][filt]['end_idx']
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
    fnames = fnames_sorted[start_idx:end_idx] #Mab is found in block 5-7 for H-band



    # find the x,y position and velocity at each time step
    pos_x, pos_y = [], []
    vec_x, vec_y = [], []
    for fname in fnames:
        # load it
        hdr = fits.open(fname)[0].header
        
        # match ephemeris time with time in fits header
        obsdate = hdr["DATE-OBS"].strip(", \n")
        obsdate = datetime.strftime(datetime.strptime(obsdate, "%Y-%m-%d"), "%Y-%b-%d")
        start_time = obsdate + " " + hdr["EXPSTART"][:5]
        ephem_line = ephem.loc[start_time]
        x_shift = float(ephem_line["sat_X"])
        y_shift = float(ephem_line["sat_Y"])
        if filt == 'H':
            exposure_starts.append(start_time)
        
        # translate from arcsec to number of pixels
        x0 = x_shift / pixscale_arcsec
        y0 = y_shift / pixscale_arcsec
        
        # find the velocity vector
        nextmin = datetime.strftime(datetime.strptime(start_time, "%Y-%b-%d %H:%M") + timedelta(minutes=2), "%Y-%b-%d %H:%M")
        ephem_nextmin = ephem.loc[nextmin]
        xi = float(ephem_nextmin["sat_X"]) / pixscale_arcsec
        yi = float(ephem_nextmin["sat_Y"]) / pixscale_arcsec
        
        dx, dy = xi - x0, yi - y0
        
        # note: x has flipped sign
        pos_x.append(-x0)
        pos_y.append(y0)
        vec_x.append(-dx)
        vec_y.append(dy)

    # print the average pixel motion per frame
    #print(vec_x[-1], vec_y[-1])
    print(f'Average motion of {code} during a single exposure (pixels): {np.mean(vec_x)}, {np.mean(vec_y)}')
    
    if filt == 'H':
        ephem_h = np.array([pos_x, pos_y])
    elif filt == 'K':
        ephem_k = np.array([pos_x, pos_y])
    else:
        raise ValueError(f'Filter {filt} not recognized')

# get the first frame and find the pixel coordinates of planet center using planetnav
hdul = fits.open(fnames[0])
header = hdul[0].header
data = hdul[0].data
obs_time = header['DATE-OBS'] + ' ' + header['EXPSTART'][:-4]
start_time = datetime.strptime(obs_time, '%Y-%m-%d %H:%M:%S')
end_time = start_time + timedelta(minutes=1)
epochs = {'start':obs_time, 'stop':end_time.strftime('%Y-%m-%d %H:%M:%S'), 'step':'1m'}

obj = Horizons(id='799', location=obscode, epochs=epochs) #Uranus, Keck
ephem = obj.ephemerides()
d_AU = ephem['delta'][0]*u.au
dist = d_AU.to(u.km).value
pixscale_km = dist*np.tan(np.deg2rad(pixscale_arcsec/3600.))
#print(f'pixel scale is {pixscale_km} km')

ura = Body('Uranus', epoch=obs_time, location='568') #Keck
nav = Nav(data, ura, pixscale_arcsec)
#ephem = nav.ephem
dx_canny, dy_canny, _, _ = nav.colocate(mode='canny',
                    low_thresh=1e-5, 
                    high_thresh=0.01, 
                    sigma=5, 
                    tb=1600, 
                    a=0.1, 
                    beam=0.1,
                    diagnostic_plot=False)

ctr = (data.shape[0]/2 + dx_canny, data.shape[1]/2 + dy_canny)


fs = 14
origins = np.array([pos_x, pos_y])
fig, ax = plt.subplots(1,1, figsize = (5.5,8))

ax.imshow(data, origin='lower', vmin=0, vmax = 500, cmap='Greys_r')

colors = ['red', 'blue']
for i, positions in enumerate([ephem_k, ephem_h]):
    pos_x = np.array(positions[0]) + ctr[0]
    pos_y = np.array(positions[1]) + ctr[1]
    print(np.nonzero(pos_y > 0))
    ax.scatter(pos_x, pos_y, color = colors[i], label = f'{filters[i]}-band')
#ax.quiver(*origins, vec_x, vec_y)
    if filters[i] == 'H':
        ax.text(pos_x[0], pos_y[0], exposure_starts[0], color = 'white', ha='left', va='bottom')
        ax.text(pos_x[-1], pos_y[-1], exposure_starts[-1], color = 'white', ha='left', va='bottom')

ax.set_ylim([0, 600])
ax.set_xlim([0, 400])
ax.set_xlabel('Pixels', fontsize = fs)
ax.set_ylabel('Pixels', fontsize = fs)
ax.tick_params(which='both', labelsize = fs - 2)
ax.legend()

# add annotation for Miranda and Puck
ax.annotate('Puck', (233, 275), xytext=(260,230), 
            color='cyan', 
            fontsize=fs, 
            arrowprops={'width':1, 'color':'cyan', 'headwidth':6, 'headlength':10}, 
            transform=ax.transAxes)
ax.annotate('Miranda', (75, 93), xytext=(90,40), 
            color='cyan', 
            fontsize=fs, 
            arrowprops={'width':1, 'color':'cyan', 'headwidth':6, 'headlength':10}, 
            transform=ax.transAxes)

plt.tight_layout()
fig.savefig(paths.figures / f"motion_on_detector_{code}.png", dpi=300)
#plt.show()
plt.close()

