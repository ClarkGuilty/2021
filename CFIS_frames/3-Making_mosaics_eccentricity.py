#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 19:27:59 2021

@author: Javier Alejandro Acevedo Barroso
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import wcs
from astropy.io import fits
from astropy.visualization import simple_norm


from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u

from numpy.random import Generator, PCG64
rng = Generator(PCG64(3))

#%% #Taken from https://github.com/esavary/Visualisation-tool/blob/master/visualisation_1band.py

def background_rms_image(cb, image):
    xg, yg = np.shape(image)
    cut0 = image[0:cb, 0:cb]
    cut1 = image[xg - cb:xg, 0:cb]
    cut2 = image[0:cb, yg - cb:yg]
    cut3 = image[xg - cb:xg, yg - cb:yg]
    l = [cut0, cut1, cut2, cut3]
    m = np.mean(np.mean(l, axis=1), axis=1)
    ml = min(m)
    mm = max(m)
    if mm > 5 * ml:
        s = np.sort(l, axis=0)
        nl = s[:-1]
        std = np.std(nl)
    else:
        std = np.std([cut0, cut1, cut2, cut3])
    return std

def scale_val(image_array):
    if len(np.shape(image_array)) == 2:
        image_array = [image_array]
    vmin = np.min([background_rms_image(8, image_array[i]) for i in range(len(image_array))])
    xl, yl = np.shape(image_array[0])
    box_size = 18  # in pixel
    xmin = int((xl) / 2. - (box_size / 2.))
    xmax = int((xl) / 2. + (box_size / 2.))
    vmax = np.max([image_array[i][xmin:xmax, xmin:xmax] for i in range(len(image_array))])
    # return np.min([vmin,vmax]), np.max([vmin,vmax])
    return vmin, vmax

#%%

frame_name = "CFIS.012.243.r"
cat_name = "cleaned.csv"
fontsize = 26
list_of_frames_file='list_of_frames'
# list_of_frames_file='list_of_objects.test'


#%%
bin_ranges = np.arange(0.6,0.8,0.05,dtype=np.float16)
size_of_mosaic = 8    
delta = bin_ranges[1] - bin_ranges[0]
with open(list_of_frames_file, 'r') as file:
    for file in file:
        if file[-2] == 'u':
            continue
        frame_name=file[:-1]
        print(frame_name)
        with fits.open(frame_name+'.fits') as hdul:
            w = wcs.WCS(hdul[0].header)
            df = pd.read_csv(frame_name+'/'+cat_name)
            for low_mag in bin_ranges:
                sel_df = df[np.logical_and(df['ELLIPTICITY'] > low_mag, df['ELLIPTICITY'] < low_mag+delta)]
                if np.sqrt(len(sel_df)) < 2:
                    continue
                size_of_mosaic = int(np.sqrt(len(sel_df)))
                iis = rng.choice(range(len(sel_df)), (size_of_mosaic,size_of_mosaic),replace=False)
                np.savetxt('random_numbers/'+file+'_ellip.csv',iis,delimiter=',')
                fig = plt.figure(figsize=(8,8), dpi=400)
                gridspec_kw = {'wspace':0.0001,
                               'hspace':0.01}
                
                ax_array = fig.subplots(size_of_mosaic, size_of_mosaic,
                                        squeeze=True,gridspec_kw=gridspec_kw)
                
                fig.suptitle(frame_name+'          '+str(low_mag)+' < MAG < '+ str(low_mag+delta),y=0.908,fontsize=20)
                for i in range(iis.shape[0]):
                    for j in range(iis.shape[1]):
                        # print(i,j)
                        ra, dec = sel_df.iloc[iis[i,j]][['ALPHA_J2000', 'DELTA_J2000']]
                        coords = SkyCoord(ra, dec, unit="deg",obstime=hdul[0].header['DATE'])
                        cutout = Cutout2D(hdul[0].data, coords, (64, 64), wcs=w, copy=True)
                        image = cutout.data
                        ax_array[i,j].tick_params(left=False,
                                bottom=False,
                                labelleft=False,
                                labelbottom=False)
                        
                        #Preparing the scaling. Based on Visualisation-tool by esavary.
                        scale_min, scale_max = scale_val(image)                
                        factor = np.log10(scale_max - scale_min)
                        indices0 = np.where(image < scale_min)
                        indices1 = np.where((image >= scale_min) & (image <= scale_max))
                        indices2 = np.where(image > scale_max)
                        image[indices0] = 0.0
                        image[indices2] = 1.0
                    
                        image[indices1] = np.log10(image[indices1]) / (factor * 1.0)
                        ax_array[i,j].imshow(image, cmap='gray', origin='lower')
                fig.savefig("mosaics/"+frame_name+'_'+str(low_mag)+'-'+str(low_mag+delta)+'.png', dpi=400,bbox_inches='tight',pad_inches = 0)
                fig.savefig(frame_name+"/ellip_"+str(low_mag)+'-'+str(low_mag+delta)+'.png',bbox_inches='tight',pad_inches = 0)
        
