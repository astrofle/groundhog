
"""
"""

import numpy as np

from groundhog import sd_fits_utils
from groundhog.scan import Scan


class SDFITS:
    
    __name__ = "SDFITS"
    
    def __init__(self, table, header):
        
        self.table = table
        self.header = header 
        #self.unique = sd_fits_utils.parse_sdfits(table)
        #self._scans = table["SCAN"]
        #self._ifnum = unq.ifnum
        #self._cal = unq.cal

    
    def get_scans(self, scans, ifnum=None, sig=None, cal=None, plnum=None, intnum=None):
        """
        Returns the rows in the SDFITS table for the requested scan numbers.
        
        Parameters
        ----------
        scans : int or array_like
            Scans to get.
        
        Returns
        -------
        scan : `groundhog.Scan`
            Rows with the selected scans.
        """
        
        mask = np.isin(self.table["SCAN"], scans)
        
        if ifnum is not None:
            mask = mask & (self.table["IFNUM"] == ifnum)
            
        if sig is not None:
            mask = mask & (self.table["SIG"] == sig)
            
        if cal is not None:
            mask = mask & (self.table["CAL"] == cal)
            
        if plnum is not None:
            mask = mask & np.isin(self.table["PLNUM"], plnum)
        
        table_scans = self.table[mask]
        
        if intnum is not None:
            if not hasattr(intnum, "__len__"):
                intnum = [intnum]
            table_scans = table_scans[intnum]
        
        scan = Scan(table_scans)
        
        return scan
    
    
    def remove_edge_chans(self, frac=0.2, chan0=None, chanf=None):
        """
        Removes the edge channels of the SDFITS DATA table.
        
        Parameters
        ----------
        frac : float, optional
            Fraction of the edge channels to remove.
            It will remove half of this value at each end of the spectra,
            i.e., if `frac=0.2` it will remove 10% of the channels on the 
            left and 10% of the channels on the right.
        """
        
        data = self.table['DATA']
        if len(data.shape) == 1:
            nchan = data.shape[0]
        elif len(data.shape) == 2:
            nchan = data.shape[1]
        if chan0 is None:
            chan0 = int(nchan*frac/2)
        if chanf is None:
            chanf = int(nchan - nchan*frac/2)
        if len(data.shape) == 1:
            data = data[chan0:chanf]
        elif len(data.shape) == 2:
            data = data[:,chan0:chanf]
            
        new_table = sd_fits_utils.update_table_column(self.table, 'DATA', data)
        self.table = new_table
        
            
