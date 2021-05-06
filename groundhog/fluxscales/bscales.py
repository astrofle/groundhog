"""
Brightness scales.
"""

import numpy as np

from astropy import units as u
from astropy import constants as ac

from groundhog import utils
from groundhog import telescopes


def jy2k(freq, eta_a_low_freq=telescopes.gbt['aperture efficiency'],
         surf_rms=telescopes.gbt['surface rms']):
    """
    Conversion factor between Jy and K.
    
    Parameters
    ----------
    freq : `~astropy.units.Quantity`
        Frequency at which to evaluate the conversion factor.
    eta_a_low_freq : float, optional
        Low frequency aperture efficiency of the telescope.
    surf_rms : `~astropy.units.Quantity`, optional
        Surface root-mean-squared error.
        
    Returns
    -------
    jy2k : `~astropy.units.Quantity`
        Conversion factor in K/Jy.
    """
    
    lmbd = (ac.c/freq)
    eta_a = utils.ruze(lmbd, eta_a_low_freq, surf_rms)
    # Specific gain: (2k/Ap)
    gain = 2.84
    
    return gain*eta_a*u.K/u.Jy
