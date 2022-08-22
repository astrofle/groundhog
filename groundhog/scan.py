"""
Scan object.
"""

import numpy as np

from groundhog import spectral_axis


class Scan():
    
    __name__ = "Scan"
    
    def __init__(self, array):
        
        self.array = array
        self.data = array["DATA"]
        self.tsys = array["TSYS"]
    
    
    def average(self):
        """
        Averages the integrations in a scan along the time axis.
        """
        
        tint = self.array["EXPOSURE"]
        dnu = self.array["CDELT1"]
        data = np.ma.masked_invalid(self.data)
        data_avg = np.ma.average(data, axis=0, weights=dnu*tint*np.power(self.tsys, -2.))
        self.data = data_avg
        #freq_avg = np.average(self.freq, axis=0, weights=dnu*tint*np.power(self.tsys, -2.))
        #self.freq = freq_avg

        # Update the array columns.
        self.array["EXPOSURE"] = tint.sum()
        #self.array["DATA"] = data_avg
        self.array["CRVAL1"] = np.ma.average(self.array["CRVAL1"], axis=0, weights=dnu*tint*np.power(self.tsys, -2.))
        self.array["CRPIX1"] = np.ma.average(self.array["CRPIX1"], axis=0, weights=dnu*tint*np.power(self.tsys, -2.))
        self.array["CDELT1"] = np.ma.average(self.array["CDELT1"], axis=0, weights=dnu*tint*np.power(self.tsys, -2.))
        self.array["VFRAME"] = np.ma.average(self.array["VFRAME"], axis=0, weights=dnu*tint*np.power(self.tsys, -2.))


    def update_xaxis(self):
        """
        """

        self.xaxis = spectral_axis.compute_spectral_axis(self.array)

        
    def get_xaxis(self):
        """
        """
        try:
            return self.xaxis
        except AttributeError:
            self.xaxis = spectral_axis.compute_spectral_axis(self.array)
            return self.xaxis
