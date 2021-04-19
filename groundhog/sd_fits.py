
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
            mask = mask & (self.table["IFNUM"] == sig)
            
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
