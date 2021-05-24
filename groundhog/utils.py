"""
Utility functions.
"""

import numpy as np

from functools import reduce


def factors(n):
    """
    Decomposes a number into its factors.
    
    Parameters
    ----------
    n : int
        Number to decompose.
    Returns
    -------
    factors : array
        Array of values into which n can be decomposed.
    """
    
    factors = set(reduce(list.__add__,
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))
    factors = np.fromiter(factors, int, len(factors))

    return factors


def get_ps_scan_pair(scan, procnum, procname):
    """
    """
    
    if np.all(procnum == 1) and procname == "OffOn":
        scan_on = scan + 1
        scan_off = scan
    elif np.all(procnum == 2) and procname == "OffOn":
        scan_on = scan
        scan_off = scan - 1
    elif np.all(procnum == 1) and procname == "OnOff":
        scan_on = scan
        scan_off = scan + 1
    elif np.all(procnum == 2) and procname == "OnOff":
        scan_on = scan - 1
        scan_off = scan
        
    return scan_on, scan_off
    

def ruze(lmbd, g0, surf_rms):
    """
    Ruze equation.
    
    Parameters
    ----------
    lmbd : float or `~astropy.units.quantity.Quantity`
        Wavelength.
    g0 : float
        Aperture efficiency at long wavelengths.
    surf_rms : float or `~astropy.units.quantity.Quantity`
        Surface rms.
    
    Returns
    -------
    eta_a : float
        Aperture efficiency at lmbd.
    """
    
    return g0*np.exp(-1.*np.power(4.*np.pi*surf_rms/lmbd, 2.))
