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

# Make test hdf5 file

Npix = 32

NCR_in = 1e-5 # CR volume density in cm^-3
B_in = 10 # B in Gauss



f2 = h5py.File('cut_fake_%d.hdf5'%(Npix),'w')

ncr = np.ones((Npix,Npix,Npix))*NCR_in
# add a line of zero ncr to be able to discern axes
ncr[:,0,:] = np.zeros_like(ncr[:,0,:])
ncr[:,:,10] = np.zeros_like(ncr[:,:,10])

Bz = np.zeros_like(ncr)+B_in/10.
Bx = np.zeros_like(ncr)
By = np.ones((Npix,Npix,Npix))*B_in
nth = np.zeros_like(ncr)+1e-6
dens = np.ones_like(ncr)

f2.create_dataset('n_CRE',data=ncr)
f2.attrs.create("%s_units" % 'n_CRE', '1/cm**3')

f2.create_dataset('n_th',data=nth)
f2.attrs.create("%s_units" % 'n_th', '1/cm**3')

f2.create_dataset('H_nuclei_density',data=dens)
f2.attrs.create("%s_units" % 'dens', '1/cm**3')

f2.create_dataset('magnetic_field_z',data=Bz)
f2.attrs.create("%s_units" % 'magnetic_field_z', 'gauss')

f2.create_dataset('magnetic_field_x',data=Bx)
f2.attrs.create("%s_units" % 'magnetic_field_x', 'gauss')

f2.create_dataset('magnetic_field_y',data=By)
f2.attrs.create("%s_units" % 'magnetic_field_y', 'gauss')

f2.close()  

#
