#!/usr/bin/env python

'''
Import the results images and do photometry on the moons using photutils
we take the standard star flux conversion factor from Paradis+23
we take the PSF wing correction from psf_wing_correction.py
we compute the error on the photometry by using many inner and outer radii


TO DO:
- add data for all the moons we want to include
- make this give ``final answers" and compile them into tables
'''

from photutils import aperture
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy import table
import paths
from scipy.interpolate import interp1d

code = 'Mab'
band = 'H'

photometry_constants_dict = {
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


stem = photometry_constants_dict[code][band]['stem']
n_coadds_raw = photometry_constants_dict[code][band]['coadds'] #apparently the "raw" images were already divided by the number of coadds so this is needed
xpos = photometry_constants_dict[code][band]['xpos'] 
ypos = photometry_constants_dict[code][band]['ypos'] 
# Paradis et al conversion factors for 2019-10-28
C1 = photometry_constants_dict[code][band]['C1']  # from cts s-1 to erg s-1 cm-2 um-1; see Paradis+23 Section 2.5.1

fname = paths.data / f'results/{stem}_{code}_2019-10-28.fits'
boxr = 50
rinner = 7
rmid = 12
router = 45
CU, CR = 1,1 # see Paradis+23 Section 2.5. this could be looked at later but for now it's ok



hdul = fits.open(fname)
itime = hdul[0].header['ITIME']/n_coadds_raw #seconds
data = hdul[0].data

# find the wing correction by importing tables from psf_wing_correction.py
wing_corrs = {}
for i in [0,1,2]:
    wing_corr_table = ascii.read(paths.data / f"tables/wing_correction_{band.upper()}_{i}.csv")
    wing_corr_interp = interp1d(wing_corr_table['Radius'], wing_corr_table['Correction'])
    wing_corrs[f"{i}"] = wing_corr_interp

def get_wing_corr(r, wing_corrs):
    '''wing_corrs is a dictionary defined above'''
    corr0 = wing_corrs["0"](r)
    corr1 = wing_corrs["1"](r)
    corr2 = wing_corrs["2"](r)
    return np.mean(np.array([corr0, corr1, corr2]))

# compute the photometry for a single inner and outer radius
ap = aperture.CircularAperture((xpos,ypos), r=rinner)
ann = aperture.CircularAnnulus((xpos,ypos), r_in=rmid, r_out=router)
apstats = aperture.ApertureStats(data, ap)
annstats = aperture.ApertureStats(data, ann)
wing_correction = get_wing_corr(rinner, wing_corrs)
bkgd = annstats.median
flux = apstats.sum - (bkgd * ap.area)
#print(f'Un-corrected total flux = {flux}')
corrected_flux = (C1/itime) * (flux / wing_correction)
#print(f'Corrected total flux = {corrected_flux} erg s-1 cm-2 um-1')
rms = annstats.std * np.sqrt(ap.area)
print(f'Approximate SNR = {flux / rms}')


# plot
doplot = True
if doplot:
    ap_patches = ap.plot(color='blue', lw=3, linestyle = '-',
                               label='Photometry aperture')
    ann_patches = ann.plot(color='red', lw=3, linestyle = '-.',
                                        label='Background annulus')
    plt.imshow(data - bkgd, origin = 'lower', vmin=-20, vmax=50, cmap="Greys_r")
    plt.xlim([xpos-boxr,xpos+boxr])
    plt.ylim([ypos-boxr,ypos+boxr])
    plt.xlabel('Pixels')
    plt.ylabel('Pixels')
    plt.savefig(paths.figures / f"example_aperture_{code}_{stem}.png", dpi=300)
    plt.show()
    plt.close()


## compute the photometry for many inner and outer radii
nsamples = 50
ri0 = 5
ri1 = 12
ro0 = ri1 + 5
ro1 = 50

# first get 100 backgrounds, 100 sums/areas, and the corresponding wing corrections
ris = np.linspace(ri0, ri1, nsamples)
ros = np.linspace(ro0, ro1, nsamples)
aps = np.array([aperture.CircularAperture((xpos,ypos), r=ri) for ri in ris])
areas = np.array([ap.area for ap in aps])
apsums = np.array([aperture.ApertureStats(data, ap).sum for ap in aps])
anns = np.array([aperture.CircularAnnulus((xpos,ypos), r_in=ri1, r_out=ro) for ro in ros])
bkgds = np.array([aperture.ApertureStats(data, ann).median for ann in anns])
wing_corrs = np.array([get_wing_corr(ri, wing_corrs) for ri in ris])

# do photometry for all inner, outer radii with the approprate wing corrections
results_table = np.empty((len(apsums), len(bkgds)))
for i in range(len(apsums)):
    apsum = apsums[i]
    area = areas[i]
    wing_correction = wing_corrs[i]
    for j, bkgd in enumerate(bkgds):
        flux = apsum - (bkgd * area)
        corrected_flux = (C1/itime) * (flux / wing_correction)
        results_table[i,j] = corrected_flux
        
final_flux = np.mean(results_table)
final_err = np.abs(np.percentile(results_table, [16, 84]) - final_flux)
print(f'{code} flux = {final_flux} erg s-1 cm-1 um-1')
print(f'Lower, upper error = {final_err} erg s-1 cm-1 um-1')
print(f'Standard deviation = {np.std(results_table)} erg s-1 cm-1 um-1')
        
fig, (ax0) = plt.subplots(1,1, figsize = (8,6))

cim0 = ax0.contourf(ros, ris, results_table, levels=12)
cbar0 = fig.colorbar(cim0)
cbar0.ax.set_ylabel('Derived Flux')
ax0.set_ylabel('inner radii')
ax0.set_xlabel('outer radii')

#ax1.hist(results_table.flatten())
#ax1.axvline(final_flux, color = 'k', linestyle = '-')
#ax1.axvline(final_err[0], color = 'k', linestyle = '--')
#ax1.axvline(final_err[1], color = 'k', linestyle = '--')
#ax1.set_ylabel('Number of Samples')
#ax1.set_xlabel('Derived Flux')

fig.savefig(paths.figures / f'photometry_vs_region_{code}_{stem}.png', dpi=300)
plt.show()
