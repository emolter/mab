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
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


code = 'Mab'
constants_dict = {
    'Mab':{
        'H':{'stem':'urh', 'coadds':6, 'C1':8.87e-17, 'xpos':187, 'ypos':395}, #urh156
        'K':{'stem':'urk', 'coadds':4, 'C1':6.46e-17, 'xpos':188, 'ypos':307}, #urk146
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

box_x = 60
box_y = 30

# parameters of the center-surround
fwhm_px = 4.0 #approx number of pixels of Keck diffraction limit
sigma_inner = fwhm_px/2.35482004503
sigma_outer = 3*sigma_inner

# parameters of the binning
zfactor = 2

for band in ['H', 'K']:
    stem = constants_dict[code][band]['stem']
    fname_real = paths.data / f'results/{stem}_{code}_2019-10-28.fits'
    fstem_experiment = paths.data / 'perturbation_experiment'
    hdul_real = fits.open(fname_real)
    data_real = hdul_real[0].data
    ypos = constants_dict[code][band]['ypos']
    xpos = constants_dict[code][band]['xpos']
    
    data_real = data_real[ypos-box_y:ypos+box_y, xpos-box_x:xpos+box_x]
    cs_conv_real = difference_of_gaussians(data_real, sigma_inner, sigma_outer)
    print(np.sum(data_real), np.sum(cs_conv_real)) #this does not conserve flux!
    ## see workaround below: just rescale by std of data and farmes

    # rebin factor of z
    cs_rebin_real=rebin(cs_conv_real, 1/zfactor) #this step does conserve flux
    detection_real = np.max(cs_rebin_real) /(np.std(cs_rebin_real))
    
    # this scaling is inconsistent with Imke's solution for upper limits
    # need to instead convert this to the nsigma of the maximum single pixel
    # or rescale Imke's version so that we compare the per-pixel RMS 
    # with the total area of the psf
    
    detection_sizes = []
    for i in range(100):
        fname = fstem_experiment / f'{stem}_{code}_2019-10-28_{i}.fits'
        data = fits.open(fname)[0].data
        data = data[ypos-box_y:ypos+box_y, xpos-box_x:xpos+box_x]
        cs_conv = difference_of_gaussians(data, sigma_inner, sigma_outer)
        cs_rebin=rebin(cs_conv, 1/zfactor) #factor of z**2 to conserve flux
        #err = np.std(cs_rebin) #by binning, we are increasing signal by factor of z^2 so snr increases by factor of z
        detection = np.max(cs_rebin) / np.std(cs_rebin)
        detection_sizes.append(detection / detection_real)
 
    fig = plt.figure(figsize = (12, 8))
    ax3 = plt.subplot2grid((2,2),(0,0))
    ax2 = plt.subplot2grid((2,2),(1,0))
    ax0 = plt.subplot2grid((2,2),(0,1))
    ax1 = plt.subplot2grid((2,2),(1,1))
    
    im0 = ax0.imshow(cs_rebin_real, origin='lower')
    ax0.set_title('Filtered and Binned Data')
    im2 = ax2.imshow(cs_rebin, origin='lower')
    ax2.set_title('Example Test Frame, Filtered and Binned')
    
    im3 = ax3.imshow(data_real, origin='lower', vmin=-30, vmax=50)
    ax3.set_title('Real Unfiltered Data')
    
    ax1.hist(detection_sizes, bins=10)
    ax1.axvline(1.0, color='k', linestyle='--', label='Panel (b) single-pixel max')
    ax1.set_xlabel('Normalized Flux of Brightest Pixel')
    ax1.set_ylabel('Number of Test Frames')
    ax1.legend(loc = 'upper right')
    
    ims = [im0, im2, im3]
    for j, ax in enumerate([ax0, ax2, ax3]):
        ax.set_xticks([])
        ax.set_yticks([])
        
        ax_divider = make_axes_locatable(ax)
        cax = ax_divider.append_axes("bottom", size="9%", pad="2%")
        cb = fig.colorbar(ims[j], orientation='horizontal', cax=cax, label='Flux (arb. units)')
        
        
    panel_labels = ['(a)', '(b)', '(c)', '(d)']
    for i, ax in enumerate([ax3, ax0, ax2, ax1]):
        
        ax.text(0.02, 0.98, panel_labels[i], ha='left', va='top', transform=ax.transAxes, fontsize = 22)
        
    
    fig.savefig(paths.figures / f'random_stack_experiment_{code}_{band}.png', dpi=300)
    #plt.show()
    plt.close()
    
    detection_sizes = np.array(detection_sizes)
    percent_higher = int(100*np.sum(detection_sizes > 1.0) / detection_sizes.size)
    with open(paths.output / f'perturbation_percent_higher_{band}.txt', 'w') as f:
        print(f"{percent_higher}", file=f)
    print(f'{percent_higher} percent of test frames contained a larger flux spike in {band} band')