"""
"""

import numpy as np

from groundhog import spectral_axis
from astropy.nddata import NDDataArray

class Scan():
    
    __name__ = "Scan"
    
    def __init__(self, table):
        
        self.table = table
        self.data = np.ma.masked_invalid(table["DATA"])
        self.tsys = table["TSYS"]
        self.get_freq()
        #self.freq = None
    
    
    def average(self):
        """
        Averages the integrations in a scan along the time axis.
        """
        
        tint = self.table["EXPOSURE"]
        dnu = self.table["CDELT1"]
        data_avg = np.ma.average(self.data, axis=0, weights=dnu*tint*np.power(self.tsys, -2.))
        self.data = data_avg
        freq_avg = np.average(self.freq, axis=0, weights=dnu*tint*np.power(self.tsys, -2.))
        self.freq = freq_avg
        
        
    def get_freq(self):
        """
        """

        self.freq = spectral_axis.compute_freq_axis(self.table)
