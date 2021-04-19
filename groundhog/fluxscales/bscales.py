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
    """
    
    lmbd = (ac.c/freq)
    eta_a = utils.ruze(lmbd, eta_a_low_freq, surf_rms)
    # Specific gain: (2k/Ap)
    gain = 2.84
    
    return gain*eta_a*u.K/u.Jy
