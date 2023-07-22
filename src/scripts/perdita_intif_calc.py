#!/usr/bin/env python

import numpy as np

'''
very simple script to compute Perdita's integrated I/F with error based on what is in
Karkoschka 01
'''

a16 = 0.033 #average Portia group reflectivity at alpha=16 degrees
phase_corr = 2.2 #from Portia group curve in other Karkoschka01 paper

radius = 15 #km
r_err = 3 #assuming this includes the uncertainty on a16

intifs = []
for r in [radius - r_err, radius, radius+r_err]:
    
    A = np.pi*r**2
    intif = A*a16*phase_corr
    intifs.append(intif)
    
print(intifs[1], intifs[1]-intifs[0], intifs[2]-intifs[1])
