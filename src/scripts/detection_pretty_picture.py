#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from astropy.io import fits
import paths

'''
Make a nicely-stretched image of the detection
Show expected motion of Mab over a single frame as a vector nearby
Show the first and second half of the observations separately as other frames?
Or else show brighter detections

"
figure caption:
Detection of Mab using the shift-and-stack technique. The red arrow indicates the average distance
and direction of Mab's motion across the detector during a single 2-minute exposure.


The average motion of Mab over a single 2-minute exposure, projected onto the detector plane,
is $\approx$6.2 pixels in the sky-North-Northeast direction. 
This is larger than Keck's diffraction-limited resolution of 40 mas (4 pixels) at 1.6 $\mu$m,
so smearing is substantial in the N-S direction, as seen in Figure \ref{fig:detection}.
"

ready to ship!
'''

code = 'Mab'

constants_dict = {
    'Mab':{
        'H':{'stem':'urh', 'coadds':6, 'xpos':187, 'ypos':395}, #urh156
        'K':{'stem':'urk', 'coadds':4, 'xpos':188, 'ypos':307}, #urk146
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

pad_x = 50
pad_y = pad_x


fs = 14
fig, axes = plt.subplots(2,1, figsize = (5.5, 8))

for i, band in enumerate(['H', 'K']):
    
    stem = constants_dict[code][band]['stem']
    coadds = constants_dict[code][band]['coadds']
    mab_loc = (constants_dict[code][band]['xpos'], constants_dict[code][band]['ypos'])

    hdul_full = fits.open(paths.data / f'results/{stem}_{code}_2019-10-28.fits')
    itime = hdul_full[0].header['ITIME'] / coadds # the 6 coadds per frame have already been accounted for, according to Imke
    data_full = hdul_full[0].data / itime
    
    ax0 = axes[i]
    cim0 = ax0.imshow(data_full, vmin = -5e-2, vmax = 15e-2, cmap='Greys_r')
    ax0_divider = make_axes_locatable(ax0)
    cax0 = ax0_divider.append_axes("right", size="7%", pad="2%")
    cbar0 = fig.colorbar(cim0, cax=cax0)
    cbar0.set_label(r'Flux Density (cts s$^{-1}$)', fontsize = fs)

    #ax0.quiver(10 + mab_loc[0], mab_loc[1], average_motion[0], average_motion[1],
    #                    angles='xy', scale_units='xy', scale=1,
    #                    color = 'red', headlength=3, headwidth=3, headaxislength=3,
    #                    label=f'Motion of {code} during exposure')
    
    ax0.set_xlim([mab_loc[0] - pad_x, mab_loc[0] + pad_x])
    ax0.set_ylim([mab_loc[1] - pad_y, mab_loc[1] + pad_y])
    #ax0.set_xlabel('Pixels')
    #ax0.set_ylabel('Pixels')
    ax0.set_xticks([])
    ax0.set_yticks([])
    
    ax0.text(0.02, 0.98, f'{band}-band', 
                color='k', ha='left', va='top', transform=ax0.transAxes,
                fontsize = fs+4, backgroundcolor='white')
    
#ax0.legend()
plt.tight_layout()
fig.savefig(paths.figures / f'detection_images_{code}.png', bbox=None)
plt.show()

'''
# test if noise is biased in one direction or the other to make Chris happy
data_nearby = data_full[mab_loc[1] - pad_y: mab_loc[1] + pad_y, mab_loc[0] - pad_x:mab_loc[0] + pad_x]

plt.imshow(data_nearby, origin = 'lower')
plt.show()

oneway = np.mean(data_nearby, axis = 0)
otherway = np.mean(data_nearby, axis = 1)
std = np.std(data_nearby)/np.sqrt(oneway.size)

fig, ax = plt.subplots(1,1, figsize = (8,6))
ax.plot(oneway, label = 'east-west')
ax.plot(otherway, label = 'north-south')
ax.axhline(0, linestyle = '--', color = 'k')
ax.axhline(-std, linestyle = '-', color = 'k', label = 'standard deviation')
ax.axhline(std, linestyle = '-', color = 'k')
ax.set_xlabel('Pixel')
ax.set_ylabel('Average flux (cts/s)')
ax.legend()
plt.show()


'''
