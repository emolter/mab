#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy import table
from photutils import aperture
import os, glob
import paths

'''
Use standard star observations to compute the flux correction
    due to the fact that we miss the wings of the PSF when we determine the moon flux

ready to ship!
'''

if not os.path.exists(paths.data / "tables"):
    os.mkdir(paths.data / "tables")


show_figures = False
for filt in ['H', 'K']:

    nsamples = 50
    ri0 = 3 #smallest moon aperture radius
    ri1 = 20 #largest moon aperture radius
    ro0 = 100 #star aperture radius to capture all star flux; also sets inner radius of noise annulus
    ro1 = 300
    noise_ann_width = 100 #width of noise annulus starting at ro
    
    star_files = [paths.data / f"reduced/2019oct28" / s for s in os.listdir(paths.data / f"reduced/2019oct28/") if s.startswith(f'hd{filt}')]
    ris = np.linspace(ri0, ri1, nsamples)
    ros = np.linspace(ro0, ro1, nsamples)
    
    
    for k, fname in enumerate(star_files):
        
        hdul = fits.open(fname)
        #print(fname)
        data = hdul[0].data
        if filt == 'H':
            cutoff = 5800
        elif filt == 'K':
            cutoff = 2500
        data[data>cutoff] = 0 #this removes bad pixels
        star_pos = np.unravel_index(np.argmax(data), data.shape)[::-1]
        #print(star_pos, np.max(data))
        
        # flux from moon inner radius
        moon_aps = np.array([aperture.CircularAperture(star_pos, r=ri) for ri in ris])
        moon_areas = np.array([ap.area for ap in moon_aps])
        moon_apsums = np.array([aperture.ApertureStats(data, ap).sum for ap in moon_aps])
        
        # flux from star inner radius
        star_aps = [aperture.CircularAperture(star_pos, r=ro) for ro in ros]
        star_areas = [star_ap.area for star_ap in star_aps]
        star_apsums = [aperture.ApertureStats(data, star_ap).sum for star_ap in star_aps]
        if __name__ == "__main__":
            print('Total star counts = ', np.mean(star_apsums))
        
        # noise
        noise_anns = [aperture.CircularAnnulus(star_pos, r_in=ro, r_out=ro+noise_ann_width) for ro in ros]
        bkgds = [aperture.ApertureStats(data, noise_ann).median for noise_ann in noise_anns]
        
        # sanity check plot showing we are on top of the star
        moon_ap_patches = moon_aps[-1].plot(color='white', lw=2,
                                   label='Photometry aperture')
        noise_ann_patches = noise_anns[0].plot(color='red', lw=2,
                                            label='Background annulus')
        plt.imshow(data - bkgds[0], origin = 'lower')
        if show_figures:
            plt.show()
        plt.close()
        
        # compile all possible combinations of moon and star aperture sizes
        results_table = np.empty((len(moon_apsums), len(bkgds)))
        for i in range(len(moon_apsums)):
            moon_apsum = moon_apsums[i]
            moon_area = moon_areas[i]
            for j, bkgd in enumerate(bkgds):
                star_apsum = star_apsums[i]
                star_area = star_areas[i]
                star_flux = star_apsum - (bkgd * star_area)
                moon_flux = moon_apsum - (bkgd * moon_area)
                
                results_table[i,j] = moon_flux / star_flux
                
        fig, ax = plt.subplots(1,1, figsize = (8,8))
        cim = ax.contourf(ros, ris, results_table, levels=12)
        cbar = fig.colorbar(cim)
        cbar.ax.set_ylabel('Derived Flux Ratio')
        ax.set_ylabel('inner radii')
        ax.set_xlabel('outer radii')
        fig.savefig(paths.figures / f'psf_wing_vs_region_{filt}_{k}.png', dpi=300)
        if show_figures:
            plt.show()
        plt.close()
        
        # the plots show that the results are very insensitive to the choice of outer radius
        # which is a good thing. So we can compile a results table based on ro = 200
        results_final = results_table[:,int(len(bkgds)/2)]
        table_out = table.Table({'Radius':list(ris), 'Correction':list(results_final)})
        table_out.write(paths.data / f"tables/wing_correction_{filt}_{k}.csv", overwrite=True)
        
        
    