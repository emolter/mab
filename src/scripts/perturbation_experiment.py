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

TO DO:
- final decision about units and scaling of the filtered images
    - is it a good idea to try to conserve flux?
- put all the frames of the image on the same color scale, add a colorbar
- add panel labels, (a) to (d)
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pylanetary.utils import rebin
from scipy.ndimage import gaussian_filter
from scipy.signal import medfilt2d, convolve2d
import paths
from skimage.filters import difference_of_gaussians

fname_real = paths.data / 'results/urh_Mab_2019-10-28.fits'
fstem_experiment = paths.data / 'perturbation_experiment'
hdul_real = fits.open(fname_real)
data_real = hdul_real[0].data
ypos = 507 
xpos = 683
box_x = 60
box_y = 30

# parameters of the center-surround
fwhm_px = 5.0 #approx number of pixels of Keck diffraction limit
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

'''
# manual center-surround filter
sz = 5
outer = 5
cs = np.zeros((sz+2*outer,sz+2*outer)) - 0.1
cs[outer:sz+outer,outer:sz+outer] = 1.0
cs/=np.sum(cs)
#plt.imshow(cs)
#plt.show()
cs_conv_real = convolve2d(data_real, cs, mode='same')
'''

cs_conv_real = difference_of_gaussians(data_real, sigma_inner, sigma_outer)
print(np.sum(data_real)/np.sum(cs_conv_real)) #this does not conserve flux
# how to get this to conserve flux? it needs to, if we want this to stand as an SNR of our detection
print(np.std(data_real)/np.std(cs_conv_real))
# maybe we do not want this to stand in for an SNR, though - is that bad statistics?

# rebin factor of z
cs_rebin_real=zfactor**2*rebin(cs_conv_real, 1/zfactor) #factor of z**2 to conserve flux
detection_real = np.max(cs_rebin_real)/np.std(data_real)
#print(detection_real)
#plt.imshow(cs_rebin_real, origin='lower')
#plt.show()

detection_sizes = []
for i in range(100):
    fname = fstem_experiment / f'Mab_2019-10-28_{i}.fits'
    data = fits.open(fname)[0].data
    data = data[ypos-box_y:ypos+box_y, xpos-box_x:xpos+box_x]
    cs_conv = difference_of_gaussians(data, sigma_inner, sigma_outer)
    cs_rebin=zfactor**2*rebin(cs_conv, 1/zfactor) #factor of z**2 to conserve flux
    detection = np.max(cs_rebin)/np.std(data)
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
ax1.axvline(detection_real, color='k', linestyle='--', label=f'Detection of Mab; SNR ~{np.round(detection_real,1)}')
ax1.set_xlabel('SNR of Brightest Pixel')
ax1.set_ylabel('Number of Test Frames')
ax1.set_yticks(np.linspace(0,18,7))
ax1.legend()

for ax in [ax0, ax1, ax3]:
    ax.set_xticks([])
    ax.set_yticks([])

fig.savefig(paths.figures / 'random_stack_experiment.png', dpi=300)
plt.show()
