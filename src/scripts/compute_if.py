#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

'''
Start with the best photometry of a moon in erg s-1 cm-1 um-1 units,
determine the integrated I/F,
determine the I/F for some assumed radii, like 6 km and 12 km

use the Jupyter notebooks from the JWST proposal for help
'''

constants_dict = {
    'Mab':{
        'H':{'flux':5.2e-16, 'err':1.7e-16, }, 
        'K':{'flux':3.6e-16, 'err':0.9e-16, }, 
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