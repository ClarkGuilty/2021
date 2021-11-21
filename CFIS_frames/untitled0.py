#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 12:49:32 2021

@author: Javier Alejandro Acevedo Barroso
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

table = pd.read_csv('CFIS.202.260.r/output.cat', skiprows=(5), sep=" "
                    , usecols=(1,2,3,4), header=None,
        names = ['MAG_ISO', 'MAGERR_ISO', 'MAG_BEST', 'MAGERR_BEST'])
print(table.head())

table.hist('MAG_BEST',)
plt.xlim(0,-19)
#%%

from astropy.io import fits

with fits.open('CFIS.012.243.r.fits') as hdul:
    hdul[0].header





