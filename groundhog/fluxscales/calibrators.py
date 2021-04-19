"""
Coefficients to compute radio SED of calibrator sources.
"""

import numpy as np

from astropy import units as u

from groundhog.fluxscales import bscales


cal_coefs = {'Perley-Butler 2017':{'3C286': [0.03570, -0.1798, -0.4507, 1.2481],
                                   '3C348': [0.00000, -0.0951, -1.0247, 1.8298],
                                   '3C353': [-0.0732, -0.0998, -0.6938, 1.8627]}} # Perley-Butler 2017


def compute_sed(freq, scale, source, units='Jy'):
    """
    """
    
    coefs = cal_coefs[scale][source]
    nu = freq.to(u.GHz).value
    
    # Calibrator flux density in Jy.
    snu = np.power(10., np.polyval(coefs, np.log10(nu)))*u.Jy
    
    if 'K' in units:
        conv = bscales.jy2k(freq)
        snu *= conv
    
    return snu.to(units)
        
