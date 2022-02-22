"""
"""

import numpy as np

from astropy import units as u
from astropy import constants as ac


def compute_freq_axis(table, chstart=1, chstop=-1, apply_doppler=True):
    """
    """

    # Copied from GBT gridder: 
    # https://github.com/nrao/gbtgridder/blob/master/src/get_data.py
    
    #shape = table.field('data').shape
    shape = table['DATA'].shape
    if len(shape) == 1:
        ax = 0
    elif len(shape) == 2:
        ax = 1
    
    if chstop == -1:
        chstop = shape[ax] + 1
        
    freq = np.zeros(shape)

    crv1 = table['CRVAL1']*u.Hz
    cd1 = table['CDELT1']*u.Hz
    crp1 = table['CRPIX1']
    vframe = table['VFRAME']*u.m/u.s

    # Observatory redshift.
    beta = vframe/ac.c
    
    # Doppler correction.
    doppler = np.ones_like(beta)
    if apply_doppler:
        doppler = np.sqrt((1.0 + beta)/(1.0 - beta))
    
    # Full frequency axis in doppler tracked frame from first row.
    # FITS counts from 1, this indx refers to the original axis, before chan selection.
    indx = np.arange(chstop - chstart) + chstart
    if ax == 1:
        indx = np.tile(indx, (shape[0],1))
        freq[:,:] = (crv1[:,np.newaxis] + cd1[:,np.newaxis]*(indx - crp1[:,np.newaxis]))*doppler[:,np.newaxis]
    else:
        freq = (crv1 + cd1*(indx - crp1))*doppler
        
    return freq*u.Hz
