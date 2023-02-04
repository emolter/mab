#!/usr/bin/env python

'''
Compile the results of the perturbation experiment
Show that the shift-and-stack technique does not introduce spurious sources
    Show that the perturbed frames have a similar noise level, but no sources
        at the same SNR
    Automate this most likely by smoothing to the Keck beam size. 
        This inflates the single-pixel strength of Mab
        in the real frame while decreasing the strength of any noise spikes
Making the frames was done at the bottom of shift_stack_uranus.py

"The expected position of Mab was shifted by a random x,y offset 
    of mean amplitude 20 pixels in each frame, then the shift-and-stack algorithm was re-run
    assuming these new randomly-selected shifts, and the output frame saved. 
    This experiment was repeated 100 times with different random offsets applied. 
    We then searched for point sources in these 100 test frames 
    by applying a center-surround filter with an FWHM of 1.5 pixels\footnote{this corresponds
    roughly to the diffraction-limited beam size of Keck in H-band}, then binned the filtered
    images by a factor of 2 in each direction. Finally, we computed the amplitude 
    of the brightest pixel in each binned frame in units of the RMS noise of the input image; 
    this yielded a rough estimate of the signal-to-noise ratio of that peak. 
    Identical filtering, binning, and SNR calculations were run on the real data. 
    The detection of Mab in the real data was found to have a higher SNR than any peak 
    in any of the 100 test frames; see Figure \ref{fig:randomstack}."

To do:
make it do both filters in a loop in order to make both figs that go in the paper
but otherwise it's ready to ship
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pylanetary.utils import rebin
from scipy.ndimage import gaussian_filter
from scipy.signal import medfilt2d, convolve2d
import paths
from skimage.filters import difference_of_gaussians

code = 'Mab'
band = 'H'

constants_dict = {
    'Mab':{
        'H':{'stem':'urh', 'xpos':187, 'ypos':395}, #urh156
        'K':{'stem':'urk', 'xpos':188, 'ypos':307}, #urk146
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

stem = constants_dict[code][band]['stem']
fname_real = paths.data / f'results/{stem}_{code}_2019-10-28.fits'
fstem_experiment = paths.data / 'perturbation_experiment'
hdul_real = fits.open(fname_real)
data_real = hdul_real[0].data
ypos = constants_dict[code][band]['ypos']
xpos = constants_dict[code][band]['xpos']
box_x = 60
box_y = 30

# parameters of the center-surround
fwhm_px = 4.0 #approx number of pixels of Keck diffraction limit
sigma_inner = fwhm_px/2.35482004503
sigma_outer = 3*sigma_inner

# parameters of the binning
zfactor = 2

#plt.imshow(data_real, origin='lower', vmin = -10, vmax = 50)
#plt.xlim([xpos-box_x, xpos+box_x])
#plt.ylim([ypos-box_y, ypos+box_y])
#plt.show()


data_real = data_real[ypos-box_y:ypos+box_y, xpos-box_x:xpos+box_x]
#data_filtered = gaussian_filter(data_real, sigma=sigma)
#data_med = medfilt2d(data_real, kernel_size=3)


cs_conv_real = difference_of_gaussians(data_real, sigma_inner, sigma_outer)
## print(np.sum(data_real)/np.sum(cs_conv_real)) #this does not conserve flux
## how to get this to conserve flux? 
## it needs to, if we want this to stand as an SNR of our detection
## see workaround below: we can scale both the signal and the noise by the bin size
## and then we can take the noise level from the binned data

# rebin factor of z
cs_rebin_real=zfactor**2*rebin(cs_conv_real, 1/zfactor) #factor of z**2 to conserve flux
detection_real = np.max(cs_rebin_real)/(np.std(cs_rebin_real) / zfactor) #by binning, we are increasing signal by factor of z^2 so snr increases by factor of z

detection_sizes = []
for i in range(100):
    fname = fstem_experiment / f'{stem}_{code}_2019-10-28_{i}.fits'
    data = fits.open(fname)[0].data
    data = data[ypos-box_y:ypos+box_y, xpos-box_x:xpos+box_x]
    cs_conv = difference_of_gaussians(data, sigma_inner, sigma_outer)
    cs_rebin=zfactor**2*rebin(cs_conv, 1/zfactor) #factor of z**2 to conserve flux
    err = np.std(cs_rebin) / zfactor #by binning, we are increasing signal by factor of z^2 so snr increases by factor of z
    detection = np.max(cs_rebin)/err
    detection_sizes.append(detection)
 
fig = plt.figure(figsize = (12, 8))
ax3 = plt.subplot2grid((2,2),(0,0))
ax2 = plt.subplot2grid((2,2),(1,0))
ax0 = plt.subplot2grid((2,2),(0,1))
ax1 = plt.subplot2grid((2,2),(1,1))

ax0.imshow(cs_rebin_real, origin='lower')
ax0.set_title('Filtered and Binned Data')
ax2.imshow(cs_rebin, origin='lower')
ax2.set_title('Example Test Frame, Filtered and Binned')

ax3.imshow(data_real, origin='lower', vmin=-30, vmax=50)
ax3.set_title('Real Unfiltered Data')

ax1.hist(detection_sizes, bins=10)
ax1.axvline(detection_real, color='k', linestyle='--', label=f'Detection of {code}; SNR ~{np.round(detection_real,1)}')
ax1.set_xlabel('SNR of Brightest Pixel')
ax1.set_ylabel('Number of Test Frames')
#ax1.set_yticks(np.linspace(0,18,7))
ax1.legend(loc = 'upper right')

for ax in [ax0, ax2, ax3]:
    ax.set_xticks([])
    ax.set_yticks([])
    
panel_labels = ['(a)', '(b)', '(c)', '(d)']
for i, ax in enumerate([ax3, ax0, ax2, ax1]):
    
    ax.text(0.02, 0.98, panel_labels[i], ha='left', va='top', transform=ax.transAxes, fontsize = 22)
    

fig.savefig(paths.figures / f'random_stack_experiment_{code}_{band}.png', dpi=300)
plt.show()

detection_sizes = np.array(detection_sizes)
fraction_higher = np.sum(detection_sizes > detection_real) / detection_sizes.size
print(f'{100*fraction_higher} percent of test frames contained a larger flux spike')