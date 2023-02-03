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
'''

mab_loc = (683, 503)
mab_loc_firsthalf = mab_loc
#mab_loc_secondhalf = (767, 656) #these are not the correct numbers
pad_x = 50
pad_y = pad_x
average_motion = (-1.5348859518691143, 6.0000087209428985) #per exposure, from mab_path_on_detector.py

hdul_full = fits.open(paths.data / 'results/urh_Mab_2019-10-28.fits')
mab_loc_firsthalf = mab_loc
hdul_firsthalf = fits.open(paths.data / 'results/urh_Mab_2019-10-28_firsthalf.fits')
hdul_secondhalf = fits.open(paths.data / 'results/urh_Mab_2019-10-28_secondhalf.fits')
#print(hdul_secondhalf[0].header['EXPSTART'])

itime = hdul_full[0].header['ITIME'] / 6 # the 6 coadds per frame have already been accounted for, according to Imke
data_full = hdul_full[0].data / itime 
data_firsthalf = hdul_firsthalf[0].data
data_secondhalf = hdul_secondhalf[0].data


fig, ax0 = plt.subplots(1,1, figsize = (8, 8))

cim0 = ax0.imshow(data_full, vmin = -5e-2, vmax = 10e-2, cmap='Greys_r')
#ax1.imshow(data_firsthalf, vmin = -30, vmax = 50, cmap='Greys_r')
#ax2.imshow(data_secondhalf, vmin = -30, vmax = 50, cmap='Greys_r')
ax0_divider = make_axes_locatable(ax0)
cax0 = ax0_divider.append_axes("right", size="7%", pad="2%")
cbar0 = fig.colorbar(cim0, cax=cax0, label=r'Flux Density (cts s$^{-1}$)')

ax0.quiver(10 + mab_loc[0], mab_loc[1], average_motion[0], average_motion[1],
                    angles='xy', scale_units='xy', scale=1,
                    color = 'red', headlength=3, headwidth=3, headaxislength=3,
                    label='Motion of Mab during exposure')

for ax in [ax0]:
    
    ax.set_xlim([mab_loc[0] - pad_x, mab_loc[0] + pad_x])
    ax.set_ylim([mab_loc[1] - pad_y + 2, mab_loc[1] + pad_y + 2])
    #ax.set_xticks([])
    #ax.set_yticks([])
    ax.set_xlabel('Pixels')
    ax.set_ylabel('Pixels')
    
ax0.legend()
fig.savefig(paths.figures / 'detection_image.png')
plt.show()