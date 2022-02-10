#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 17:10:01 2021

@author: gin
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LogNorm
from astropy.io import fits

#%%





#%%

# Plot polaris output
savedir = '/home/gin/projects/FIRE/POLARIS/polarisfire/mock_example/results/fake_16/'
folder = savedir+'data/'

os.chdir(folder)

# structure of file is :'I Q U V null Nth I Q U V lambda^2xkappa Ncr'

hdu = fits.open(folder+'polaris_detector_nr0001.fits')
dat2 = hdu[0].data[0][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('Stokes I')
plt.savefig(savedir+'plots/I_d01.png')
plt.show()

print(dat2[0][0], 'Jy/pix')

# Angle of pixel
pix_ang = hdu[0].header['CDELT1A'] #deg
print('Pix size', pix_ang, hdu[0].header['CUNIT1A'])

hdu = fits.open(folder+'polaris_detector_nr0002.fits')
dat2 = hdu[0].data[0][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('Stokes I')
plt.savefig(savedir+'plots/I_d02.png')
plt.show()

print(dat2[0][0], 'Jy/pix')

hdu = fits.open(folder+'polaris_detector_nr0003.fits')
dat2 = hdu[0].data[0][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('Stokes I')
plt.savefig(savedir+'plots/I_d03.png')
plt.show()

print(dat2[0][0], 'Jy/pix')

hdu = fits.open(folder+'polaris_detector_nr0004.fits')
dat2 = hdu[0].data[0][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('Stokes I')
plt.savefig(savedir+'plots/I_d04.png')
plt.show()

print(dat2[0][0], 'Jy/pix')

print('-----Stokes QU------')

print('1st detector')
# QU
hdu = fits.open(folder+'polaris_detector_nr0001.fits')
dat2 = hdu[0].data[2][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
if dat2.mean()!=0:
    plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('U')
plt.savefig(savedir+'plots/U_d01.png')
plt.show()
print('U', dat2.min(), dat2.max(), 'Jy/pix')
dat2 = hdu[0].data[1][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
if dat2.mean()!=0:
    plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('Q')
plt.savefig(savedir+'plots/Q_d01.png')
plt.show()
print('Q', dat2.min(), np.max(dat2), 'Jy/pix')

print('2nd detector')
hdu = fits.open(folder+'polaris_detector_nr0002.fits')
dat2 = hdu[0].data[2][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
if dat2.mean()!=0:
    plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('U')
plt.savefig(savedir+'plots/U_d02.png')
plt.show()
print('U', dat2.min(), dat2.max(), 'Jy/pix')
dat2 = hdu[0].data[1][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
if dat2.mean()!=0:
    plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('Q')
plt.savefig(savedir+'plots/Q_d02.png')
plt.show()
print('Q', dat2.min(), np.max(dat2), 'Jy/pix')


print('3rd detector')
hdu = fits.open(folder+'polaris_detector_nr0003.fits')
dat2 = hdu[0].data[2][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
if dat2.mean()!=0:
    plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('U')
plt.savefig(savedir+'plots/U_d03.png')
plt.show()
print('U', dat2.min(), dat2.max(), 'Jy/pix')
dat2 = hdu[0].data[1][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
if dat2.mean()!=0:
    plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('Q')
plt.savefig(savedir+'plots/Q_d03.png')
plt.show()
print('Q', dat2.min(), np.max(dat2), 'Jy/pix')


print('4rd detector')
hdu = fits.open(folder+'polaris_detector_nr0004.fits')
dat2 = hdu[0].data[2][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
if dat2.mean()!=0:
    plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('U')
plt.savefig(savedir+'plots/U_d04.png')
plt.show()
print('U', dat2.min(), dat2.max(), 'Jy/pix')
dat2 = hdu[0].data[1][0]
plt.figure(figsize=(5,6))
plt.imshow(dat2,origin='lower', cmap = 'Blues')
if dat2.mean()!=0:
    plt.colorbar(fraction=0.046, pad=0.1,orientation='horizontal')
plt.title('Q')
plt.savefig(savedir+'plots/Q_d04.png')
plt.show()
print('Q', dat2.min(), np.max(dat2), 'Jy/pix')

wavelength =  hdu[0].header['HIERARCH WAVELENGTH1']
frequency = 2.99e8/wavelength/1e6#m/s/m -> MHz
print('Wavelength (m)', wavelength, 'frequency (MHz)', frequency)

#%%
# I


c = 2.99e10

k_B = 1.38e-16 #erg/K

pc_to_cm = 3.08e18 #cm

def I_to_K(I,wavelength):
    wavelength_cm = wavelength*100
    # convert I from Jy to cgs
    return I*1e-23/(2*k_B)*wavelength_cm**2
#%%
hdu = fits.open(folder+'polaris_detector_nr0001.fits')
dat2 = hdu[0].data[0][0]

print np.unique(dat2)
print np.unique(I_to_K(dat2.flatten(),wavelength)), 'K'


hdu = fits.open(folder+'polaris_detector_nr0002.fits')
dat2 = hdu[0].data[0][0]

print np.unique(dat2)
print np.unique(I_to_K(dat2.flatten(),wavelength)), 'K'

hdu = fits.open(folder+'polaris_detector_nr0003.fits')
dat2 = hdu[0].data[0][0]

print np.unique(dat2)
print np.unique(I_to_K(dat2.flatten(),wavelength)), 'K'
#%%

datNcr = hdu[0].data[11][0]
print np.unique(datNcr)
