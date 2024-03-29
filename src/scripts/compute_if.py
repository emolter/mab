#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pylanetary.utils import *
from astroquery.jplhorizons import Horizons
from astropy.io import ascii
from astropy import table
import paths

'''
Start with the best photometry of a moon in erg s-1 cm-1 um-1 units,
determine the integrated I/F,
plot this against the Showalter value for integrated I/F

description of types of albedo: https://the-moon.us/wiki/Albedo

To do: make work for the other moons, put them in the paper too
'''

caption = r'Available color information on Mab and Perdita compared with Puck and Miranda. See Figure \ref{fig:spectrum} caption for data sources. \label{tab:color}'

constants_dict = {
    'Mab':{
        'H':{'stem':'urh', 'bp_file':paths.static / 'h.csv', 'A':1}, # for integrated I/F use area of 1 km
        'K':{'stem':'urk', 'bp_file':paths.static / 'kp.csv', 'A':1}, 
        '500':[0.5, 50, 3],
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
        'H':{'stem':'urh', 'bp_file':paths.static / 'h.csv', 'A':1}, 
        'K':{'stem':'urk', 'bp_file':paths.static / 'kp.csv', 'A':1}, 
        '500':[0.5, 54, 22],
        },
}

def quadsum(err1, err2):
    return np.sqrt(err1**2 + err2**2)
    
def ratio_with_err(moonarr, i0, i1, nround = 2):
    '''
    Parameters
    ----------
    moonarr: 2-d array ([wls, fluxes, errs])
    i0, i1: indices on the inner axis to ratio, i.e., fluxes[i0]/fluxes[i1]
    
    Returns
    -------
    color
    color_err
    '''
    flux0 = moonarr[1][i0]
    flux1 = moonarr[1][i1]
    err0 = moonarr[2][i0]
    err1 = moonarr[2][i1]
    color = flux0 / flux1
    frac_err = quadsum(err0/flux0, err1/flux1)
    color_err = color*frac_err
    return np.round(color, nround), np.round(color_err, nround)
    

# add spectra of Miranda leading, Miranda trailing, and Puck from other authors
miranda_trailing = np.genfromtxt(paths.static / 'miranda_trailing.csv', delimiter=',').T #http://dx.doi.org/10.1051/0004-6361/201321988 #phase angle was around 0.2
#miranda_leading = np.genfromtxt(paths.static / 'miranda_leading.csv', delimiter=',').T
haumea = np.genfromtxt(paths.static / 'haumea.csv', delimiter=',').T #https://doi.org/10.1051/0004-6361/201526423

# photometry from other authors
gibbard_miranda = np.array([[1.25, 1.63, 2.12],
                        [241, 194, 169],
                        [17, 20, 10]]) # phase angle 1.9; rescale to 0.02 based on Karkoschka 01 using a factor of 1.5
karkoschka_miranda = np.array([[0.34, 0.55, 0.63, 0.67, 0.89, 0.91], 
                        [468, 449, 460, 443, 434, 428],
                        [21, 22, 26, 19, 28, 18]]) #phase angle 0.034
karkoschka_puck = np.array([[0.34, 0.55, 0.63, 0.67, 0.91], 
                        [109, 115, 134, 104, 106],
                        [51, 53, 66, 29, 8]]) #phase angle 0.034
paradis_puck = np.array([[1.62, 2.12],
                        [109, 129], 
                        [3, 4]]) # phase angle 0.02

# scale Miranda up from phase angle of 1.9 to phase angle of 0.02 based on curves from Karkoschka01
gibbard_miranda[1:]*=1.5

# concatenate tables
miranda = np.concatenate([karkoschka_miranda, gibbard_miranda], axis = 1)
puck = np.concatenate([karkoschka_puck, paradis_puck], axis = 1)

# what are the Kp/H ratios of each body?
color_dict = {'Miranda':{r'0.5 $\mu$m / 1.6 $\mu$m': ratio_with_err(miranda, 0, -2)[0],
                        'error':ratio_with_err(miranda, 0, -2)[1]},
                'Puck':{r'0.5 $\mu$m / 1.6 $\mu$m': ratio_with_err(puck, 0, -2)[0], 
                        'error':ratio_with_err(puck, 0, -2)[1],}
                }
#'Kp/H':[ratio_with_err(mab, -1, -2)[0], ratio_with_err(miranda, -1, -2)[0], ratio_with_err(puck, -1, -2)[0]],
#'Kp/H err':[ratio_with_err(mab, -1, -2)[1], ratio_with_err(miranda, -1, -2)[1], ratio_with_err(puck, -1, -2)[1]],
#r'0.5 $\mu$m / 1.6 $\mu$m':[ratio_with_err(mab, 0, -2)[0], , ],
#'vis/H err':[ratio_with_err(mab, 0, -2)[1], , ,
#}

for code in ['Perdita', 'Mab']:
    # get sun-target and Earth-target distance from Horizons
    obscode = '568'
    date = '2019-10-28'
    tstart = date+' 11:00' #this just needs to be approximate
    tend = date+' 12:00'
    horizons_obj = Horizons(
        id=code,
        location=obscode,
        epochs={"start": tstart, "stop": tend, "step": "1h"},
    )
    ephem = horizons_obj.ephemerides()[0]
    target_sun_dist = ephem['r']
    target_obs_dist = ephem['delta']
    area_km = constants_dict[code]['H']['A']
    target_omega = area_km*u.km**2 / (((target_obs_dist*u.au).to(u.km))**2) #setting area to 1 km2 means that I/F function will output integrated I/F
    
    wls, fluxes, errs = [], [], []
    for band in ['H', 'K']:
        # get Keck filter transmission
        stem = constants_dict[code][band]['stem']
        with open(paths.output / f'{code}_{stem}_flux.txt', 'r') as f:
            flux = float(f.readline()) * 1e-16 # erg s-1 cm-2 um-1
        with open(paths.output / f'{code}_{stem}_fluxerr.txt', 'r') as f:
            flux_err = float(f.readline()) * 1e-16 # erg s-1 cm-2 um-1
        bp_file = constants_dict[code][band]['bp_file']
        wl, trans = np.genfromtxt(bp_file, skip_header=1, delimiter=',', usecols=(0,1)).T
        wl = wl[~np.isnan(wl)][:-1] * u.micron
        trans = trans[~np.isnan(trans)][:-1] / 100.
        bp = np.array([wl, trans])
        
        # do the I/F calculation
        wl_eff, ioverf = I_over_F(flux, bp, target_sun_dist, target_omega)
        fractional_err = flux_err / flux
        ioverf_err = ioverf*fractional_err
        
        with open(paths.output / f"{code}_{stem}_intif.txt", "w") as f:
            print(f"{int(ioverf)}", file=f)
        with open(paths.output / f"{code}_{stem}_intiferr.txt", "w") as f:
            print(f"{int(ioverf_err)}", file=f)
        
        print(f'{code} {band}-band integrated I/F = {ioverf} +/- {ioverf_err} km2 at {wl_eff} um')
        
        
        wls.append(wl_eff)
        fluxes.append(ioverf)
        errs.append(ioverf_err)

    #put Showalter point in 
    mab = np.array([wls, fluxes, errs])
    mab = np.insert(mab, 0, constants_dict[code]['500'], axis=1)
    

    color_dict[code] = {r'0.5 $\mu$m / 1.6 $\mu$m':ratio_with_err(mab, 0, -2)[0],
                        'error':ratio_with_err(mab, 0, -2)[1]}

# output the color table
bodies = list(color_dict.keys())
vis_over_h = []
vis_over_h_err = []
for body in bodies:
    vis_over_h.append(color_dict[body][r'0.5 $\mu$m / 1.6 $\mu$m'])
    vis_over_h_err.append(color_dict[body]['error'])

color_table = {'Moon':bodies, r'0.5 $\mu$m / 1.6 $\mu$m':vis_over_h, r'1$\sigma$ error':vis_over_h_err}
color_table = table.Table(color_table)
print(color_table)
ascii.write(color_table, paths.output / 'color_table.tex', 
            Writer=ascii.Latex, caption=caption, overwrite=True,
            latexdict=ascii.latex.latexdicts['AA'])


# two possible radius determinations of Mab
r1 = 12 #km
r2 = 6 #km

# plot the spectrum
fs = 14   
fig, ax = plt.subplots(1,1, figsize = (8,6))

colors = ['k', 'gray']
for i,r in enumerate([r1, r2]):
    area = np.pi*r**2
    ioverf = mab[1]/area
    ioverf_err = mab[2]/area
    
    ax.errorbar(mab[0], ioverf, yerr=ioverf_err, uplims = [False, False, True],
                linestyle = '', marker = 'o', label = f'Mab r={r} km', color = colors[i])
    print(f'Given r = {r}, H-band albedo is {mab[1][-2]/area} +/- {mab[2][-2]/area}')
#ax2.scatter(miranda_trailing[0], miranda_trailing[1], color = 'blue', label = 'Miranda Trailing')
#ax2.scatter(haumea[0], haumea[1], color = 'red', label = 'Haumea')

ax.errorbar(miranda[0], miranda[1]*1e-3, yerr = miranda[2]*1e-3,
                   color = 'blue', marker = 's', label = 'Miranda')
ax.errorbar(puck[0], puck[1]*1e-3, yerr = puck[2]*1e-3,
                   color = 'red', marker = 'D', label = 'Puck')

ax.set_xlabel(r'Wavelength ($\mu$m)', fontsize = fs)
ax.set_ylabel('Geometric Albedo', fontsize = fs)
ax.set_ylim([0.0, 0.5])

ax.legend()
plt.tight_layout()
fig.savefig(paths.figures / 'reflectance_spectrum.png', dpi=300)
#plt.show()
plt.close()
