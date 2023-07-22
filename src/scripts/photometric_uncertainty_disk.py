#!/usr/bin/env python

'''
Use the frame-by-frame flux variations in a cloudless region on the planet Uranus 
in both H and Kp band to derive uncertainty in the photometry

The night was very photometrically stable. We derived the photometric uncertainty 
by looking at the frame-by-frame surface brightness variations of cloud-free regions
of the planet. We determined the surface brightness of the planet in two 20x20 pixel regions
of the disk in both H- and Kp-band, one just north of Uranus's equator, and one just south 
of the pole. The frame-by-frame variance displayed no clear temporal trends in either H- or Kp-band;
the fluctuations looked like Gaussian noise with standard deviations of 2.1\% in H-band and 5.0\% 
in Kp-band. These numbers agree with the spread in the fluxes derived from the three standard star
frames.
'''

from shift_stack_moons.image import *
from shift_stack_moons.shift_stack_moons import *
import warnings, paths, os


constants_dict = {
    'Mab':{
        'H':{'stem':'urh', 'pos1':(550, 500), 'pos2':(735, 485), 'edge_detect':False},
        'Kp':{'stem':'urk', 'pos1':(615, 500), 'pos2':(770, 470),'edge_detect':True},
        },
    'Ophelia':{
        'H':{},
        'Kp':{},
        },
    'Cupid':{
        'H':{},
        'Kp':{},
        },
    'Perdita':{
        'H':{},
        'Kp':{},
        },
}

# find filenames
code = 'Mab'
data_dir = paths.data / 'reduced/2019oct28'
data_files = os.listdir(data_dir)
obscode = '568'
date = '2019-10-28'
tstart = date+' 00:00'
tend = date+' 23:59'
pixscale = 0.009971
boxsize = 10


for band in ['H', 'Kp']:
    print(f"Starting {code} {band}-band")
    stem = constants_dict[code][band]['stem']
    edge_detect = constants_dict[code][band]['edge_detect']
    
    fname_list = [data_dir / s for s in data_files if s.startswith(stem)]
    
    # first align them as we do for the shift-and-stack
    frames = [Image(fname).data for fname in fname_list]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if edge_detect is not None:
            frames_centered = chisq_stack(frames, edge_detect=edge_detect)
        else:
            frames_centered = frames
    
    pos1 = constants_dict[code][band]['pos1']
    pos2 = constants_dict[code][band]['pos2']
    fluxes_1, fluxes_2 = [], []
    for j, frame in enumerate(frames_centered):
        
        x0, x1 = pos1[0]-boxsize, pos1[0]+boxsize
        y0, y1 = pos1[1]-boxsize, pos1[1]+boxsize
        box1 = frames_centered[j][y0:y1, x0:x1]
        x0, x1 = pos2[0]-boxsize, pos2[0]+boxsize
        y0, y1 = pos2[1]-boxsize, pos2[1]+boxsize
        box2 = frames_centered[j][y0:y1, x0:x1]
        
        #if j == 0:
        #    plt.imshow(box1, origin = 'lower') #expect 1.38e3 in H-band
        #    plt.show()
        #    
        #    plt.imshow(box2, origin = 'lower') #expect 3.18e3 in H-band
        #    plt.show()
            
        fluxes_1.append(np.mean(box1))
        fluxes_2.append(np.mean(box2))
    
    
    frac1 = fluxes_1/np.mean(fluxes_1) 
    frac2 = fluxes_2/np.mean(fluxes_2) 
    std1 = np.std(frac1)
    std2 = np.std(frac2) 
    meanstd = np.mean([std1, std2])
    percenterr = meanstd*100
    with open(paths.output / f"{stem}_photometric_uncertainty_percent.txt", "w") as f:
        print(f"{percenterr:.1f}", file=f)
    
    fig, ax = plt.subplots(1,1, figsize = (8,6))
    ax.plot(frac1, label = 'near equator')
    ax.plot(frac2, label = 'near pole')
    ax.axhline(1 - meanstd, color = 'k', linestyle = '--', label = f'1 sigma = {percenterr:.1f}%')
    ax.axhline(1 + meanstd, color = 'k', linestyle = '--')
    ax.axhline(1, color = 'k', linestyle = '-')
    ax.set_xlabel('Frame ID')
    ax.set_ylabel('Fractional Flux Difference')
    ax.set_title(f'{band}-band variation in disk flux')
    ax.legend(loc='upper right')
    fig.savefig(paths.figures / f'disk_flux_variation_{band}.png', dpi=300)
    #plt.show()
    plt.close()
    
