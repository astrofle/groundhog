"""
Constants for different telescopes.
"""

from astropy import units as u


gbt = {'name': 'Green Bank Telescope',
       'diameter': 100*u.m,
       'surface rms': 350*u.um,
       'aperture efficiency': 0.71}


class Telescope():
    
    def __init__(self):
        pass
