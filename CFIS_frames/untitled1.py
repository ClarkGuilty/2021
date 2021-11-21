#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 01:24:34 2021

@author: Javier Alejandro Acevedo Barroso
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from astropy.io import fits
from astropy.table import Table
from scipy.stats import norm

frame_name = "CFIS.012.243.r"
cat_name = "output.fits"
# with fits.open(frame_name+'/'+cat_name) as hdul:
#     hdul[0].header
#     print(hdul['SCI'].data)
fontsize = 26

with open('list_of_objects.test', 'r') as file:
    for file in file:
        frame_name=file[:-1]
        print(frame_name)
        
        #Read table and focus on physical values.
        table = Table.read(frame_name+'/'+cat_name,
                   format='fits',)
        df_0 = table.to_pandas()
        df = df_0[np.logical_and(
    np.logical_and(df_0['MAG_AUTO'] < 30, df_0['MAG_AUTO'] > 20),
    np.logical_and(df_0['FLUX_RADIUS'] > 0, df_0['FLUX_RADIUS'] < 10))]        # plt.scatter(df.FLUX_RADIUS , df.MAG_AUTO)
    
        #Finding tentative FWHM using a histogram.
        #Mode for the location estimator, std for the dispersion.
        plt.figure()
        y, x, _ = plt.hist(df['FLUX_RADIUS'], bins=50,density=True)
        plt.xlim(0,18)
        plt.title("FWHM estimation for "+frame_name,fontsize=18)
        tentative_fwhm = x[np.where(y == y.max())][0]
        (mu, sigma) = norm.fit(df['FLUX_RADIUS'])
        plt.plot(x, norm.pdf(x, tentative_fwhm, sigma))
        print('Tentative fwhm: ', tentative_fwhm)
        print('mu: ', mu, '\t sigma:', sigma)

        #Making the preselection.
        #Little plot of the data. The redder, the more likely it is an extended source, the greener the more likely it is a psf.
        #Plot the colored data + rectangle encapsulating the preselection.
        ax = df.plot.scatter(x = 'FLUX_RADIUS', y ="MAG_AUTO",
                c = 'CLASS_STAR', s=0.8,
                figsize=(20,15), fontsize=fontsize,
                xlim=(0,8), ylim=(26.1,19.95), cmap = 'PiYG')
        ax.xaxis.get_label().set_visible(True)
        ax.tick_params(axis='x', bottom=True, labelbottom=True,labelsize=fontsize)
        ax.set_xlabel('Flux radius (half-flux)', fontsize=fontsize)
        ax.set_ylabel('MAG Aper (Kron)', fontsize=fontsize)
        
        
        mag_high = 24; mag_low = 20
        radius_low = tentative_fwhm + sigma; radius_high = 6
        greater_than = df[np.logical_and(
            np.logical_and(df['MAG_AUTO'] < mag_high, df['MAG_AUTO'] > mag_low),
            np.logical_and(df['FLUX_RADIUS'] > radius_low, df['FLUX_RADIUS'] < radius_high))]
        
        rectangle = plt.Rectangle((radius_low,mag_low), radius_high-radius_low, mag_high-mag_low, fc='none',ec="red")
        plt.gca().add_patch(rectangle)
        ax.set_title("MAG vs R for : "+frame_name+".  R = ("+str(round(radius_low,2))+', '+str(round(radius_high,2))+')', fontsize=28)

        