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
from scipy.stats import median_abs_deviation

frame_name = "CFIS.012.243.r"
cat_name = "output.fits"
# with fits.open(frame_name+'/'+cat_name) as hdul:
#     hdul[0].header
#     print(hdul['SCI'].data)
fontsize = 26
list_of_frames_file='list_of_frames'
# list_of_frames_file='list_of_objects.test'

def mad(x):
    return np.mean(np.abs(np.median(x) - x))
#%%

with open(list_of_frames_file, 'r') as file:
    totals = []
    for file in file:
        if file[-2] == 'u':
            continue
        frame_name=file[:-1]
        print(frame_name)
        
        #Read table and focus on physical values.
        table = Table.read(frame_name+'/'+cat_name,
                   format='fits',)
        df_0 = table.to_pandas()
        df = df_0[np.logical_and(
            np.logical_and(df_0['FLUX_RADIUS'] > 0, df_0['FLUX_RADIUS'] < 10),
                np.logical_and(df_0['MAG_AUTO'] < 30, df_0['MAG_AUTO'] > 19))]        # plt.scatter(df.FLUX_RADIUS , df.MAG_AUTO)
    # np.logical_and(df_0['FLUX_RADIUS'] > 0, df_0['FLUX_RADIUS'] < 1000))]        # plt.scatter(df.FLUX_RADIUS , df.MAG_AUTO)
    
    
        #Finding tentative FWHM using a histogram.
        #Mode for the location estimator, std for the dispersion.
        plt.figure()
        print(len(df['FLUX_RADIUS']))
        y, x, _ = plt.hist(df['FLUX_RADIUS'], bins=50,density=False)
        # plt.xlim(0,18)
        plt.title("FWHM estimation for "+frame_name,fontsize=18)
        tentative_fwhm = x[np.where(y == y.max())][0]
        # tentative_fwhm = np.median(df['FLUX_RADIUS'])
        # (mu, sigma) = norm.fit(df['FLUX_RADIUS'])
        mu = np.median(df['FLUX_RADIUS'])
        sigma = mad(df['FLUX_RADIUS'])
        plt.plot(x, norm.pdf(x, tentative_fwhm, sigma)*y.max())
        plt.xlabel("FLUX_RADIUS")
        plt.savefig('histograms3/'+frame_name+"_FLUX_RADIUS_hist.png",dpi = 400)
        plt.savefig(frame_name+"/fig_fwhm_estimation.png",dpi = 400)

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
        ax.set_ylabel('MAG_AUTO', fontsize=fontsize)
        
        
        mag_high = 22.5; mag_low = 20
        # ellip_high = 0.8; ellip_low = 0.6
        ellip_high = 1; ellip_low = 0
        radius_low = tentative_fwhm + sigma; radius_high = 5000
        selection_df = df_0[np.logical_and(
            np.logical_and(df_0['ELLIPTICITY'] < ellip_high, df_0['ELLIPTICITY'] > ellip_low),
            np.logical_and(
            np.logical_and(df_0['MAG_AUTO'] < mag_high, df_0['MAG_AUTO'] > mag_low),
            np.logical_and(df_0['FLUX_RADIUS'] > radius_low, df_0['FLUX_RADIUS'] < radius_high)))]
        selection_df.to_csv(frame_name+'/cleaned.csv',index=False)
        
        #plot of square to verify rectangle.
        rectangle = plt.Rectangle((radius_low,mag_low), radius_high-radius_low, mag_high-mag_low, fc='none',ec="red")
        plt.gca().add_patch(rectangle)
        ax.set_title("MAG_AUTO vs R for : "+frame_name+".  R = ("+str(round(radius_low,2))+', '+str(round(radius_high,2))+')', fontsize=28)
        
        totals.append(len(selection_df))
        print(str(len(selection_df))+str(" selected files"))
        print("## -- -- -- -- -- -- ##\n")
        plt.savefig(frame_name+"/fig_selection.png",dpi = 400)
        
        #histogram of the magnitudes of the final selection.
        selection_df = df[np.logical_and(df['FLUX_RADIUS'] > radius_low, df['FLUX_RADIUS'] < radius_high)]
        plt.figure()
        y, x, _ = plt.hist(selection_df['MAG_AUTO'], bins=int(round(len(selection_df)/400)),density=False, cumulative=True)
        # plt.xlim(0,18)
        y, x, _ = plt.hist(df['MAG_AUTO'], bins=int(round(len(selection_df)/400)),density=False, cumulative=True)
        plt.title("MAG_AUTO hist for "+frame_name,fontsize=18)
        tentative_fwhm = x[np.where(y == y.max())][0]
        mu = np.median(df['MAG_AUTO'])
        sigma = median_abs_deviation(df['MAG_AUTO'])
        plt.xlabel("MAG_AUTO")
        plt.savefig('histograms/'+frame_name+"_MAG_hist.png",dpi = 400)

        selection_df = df[np.logical_and(df['FLUX_RADIUS'] > radius_low, df['FLUX_RADIUS'] < radius_high)]
        plt.figure()
        y, x, _ = plt.hist(1-selection_df['ELLIPTICITY'], bins=int(round(len(selection_df)/400)),density=False, cumulative=True)
        # plt.xlim(0,18)        
        plt.title("ELLIPTICITY hist for "+frame_name,fontsize=18)
        tentative_fwhm = x[np.where(y == y.max())][0]

        plt.xlabel("B/A  (1 - ELLIPTICITY)")
        plt.savefig('histograms2/'+frame_name+"_ELLIPTICITY_hist.png",dpi = 400)

print('The average number of selected galaxies per frame was:', sum(totals)/len(totals))
        
        
        
        
        