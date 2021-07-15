
"""
"""

import numpy as np

from groundhog.scan import Scan
from groundhog import datared
from groundhog import sd_fits_utils


class SDFITS:
    
    __name__ = "SDFITS"
    
    def __init__(self, hdu=None, index=None):
        
        self.hdu = hdu
        self.index = index
        
        
    def load(self, filename):
        """
        """
        
        self.hdu = fits.open(filename, memmap=True)
        self.index = sd_fits_utils.build_index(hdu, ext='SINGLE DISH')
        
    
    def get_scans(self, scans, ifnum=None, sig=None, cal=None, plnum=None, fdnum=None, intnum=None):
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
        
        table = None
        
        # Select the table from the possible tables.
        if self.numtab == 1:
            table = self.table[0]
        elif self.numtab > 1:
            for i in range(self.numtab):
                if np.isin(self.table[i]["SCAN"], scans).sum() > 0:
                    if table is not None:
                        print("Scans span multiple configurations.")
                        print("This is not supported.")
                        return
                    table = self.table[i]
        
        mask = sd_fits_utils.get_table_mask(table, scans=scans, 
                                            ifnum=ifnum, sig=sig, 
                                            cal=cal, plnum=plnum)
        
        table_scans = table[mask]
        
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
        
        for i,table in enumerate(self.table):
            data = table['DATA']
            
            if len(data.shape) == 1:
                nchan = data.shape[0]
            elif len(data.shape) == 2:
                nchan = data.shape[1]
                
            if chan0 is None:
                chan0 = int(nchan*frac/2)
            if chanf is None:
                chanf = int(nchan - nchan*frac/2)
            
            # Select only inner channels.
            if len(data.shape) == 1:
                data = data[chan0:chanf]
            elif len(data.shape) == 2:
                data = data[:,chan0:chanf]
            
            # Update DATA column.
            new_table = sd_fits_utils.update_table_column(table, 'DATA', data)
            self.table[i] = new_table
            
            # Update frequency axis header.
            crp1 = table.field('crpix1')
            crp1 -= chan0
            new_table = sd_fits_utils.update_table_column(table, 'CRPIX1', crp1)
            self.table[i] = new_table
            
    
    def update_tcal(self, scan, ifnum=None, plnum=None, update_scans=None,
                    scale="Perley-Butler 2017", units="K", avgf_min=16):
        """
        Updates the TCAL column on the SDFITS table.
        `scan` must point to one of the scans of a flux density calibrator.
        
        Parameters
        ----------
        scan : int
            Scan with observations of a calibrator source.
            This will be used to derive the temperature of the
            noise diode.
        ifnum : list, optional
            Spectral windows to process.
            Will process all spectral windows by default.
        plnum : list, optional
            Polarizations to process.
            Will process all polarizations by default.
        update_scans : list, optional
            List of scans to update with the new TCAL values.
        scale : str, optional
        
        units : {'K', 'Jy'}, optional
            TCAL units.
        avgf_min : int
            Minimum number of channels to average together when
            computing the kappa factor (Eq. (14) in Winkel et al. 2012).
        """
        
        if ifnum is None:
            ifnum = np.unique(self.table['IFNUM'])
        else:
            if not hasattr(ifnum, "__len__"):
                ifnum = [ifnum]
        
        if plnum is None:
            plnum = np.unique(self.table['PLNUM'])
        else:
            if not hasattr(plnum, "__len__"):
                plnum = [plnum]
        
        # How many channels?
        shape = self.table['DATA'].shape
        if len(shape) == 1:
            ax = 0
        elif len(shape) == 2:
            ax = 1
        nchan = shape[ax]
        
        if self.table['TCAL'].shape != shape:
            # Expand the TCAL column to accomodate vectors.
            tcal_col = self.table['TCAL']
            tcal = np.tile(tcal_col[:,np.newaxis], (1,nchan))
            self.update_table_col('TCAL', tcal)
        
        # Loop over spectral windows and polarizations
        # updating TCAL.
        for ifnum_ in ifnum:
            for plnum_ in plnum:
        
                tcal = datared.get_tcal(self, scan, ifnum=ifnum, plnum=plnum, 
                                        scale=scale, units=units, avgf_min=avgf_min)
                
                mask = sd_fits_utils.get_table_mask(self.table, scans=update_scans, 
                                                    ifnum=ifnum, plnum=plnum)
                
                tcal_col = self.table['TCAL']
                tcal_col[mask] = np.tile(tcal.value, (mask.sum(),1))
                self.update_table_col('TCAL', tcal_col)
        
    
    def update_table_col(self, column_name, column_vals, tablenum=None):
        """
        Updates `column_name` with `column_vals` in the SDFITS table.
        
        Parameters
        ----------
        column_name : str
            Name of the column to update.
        column_vals : np.ndarray
            Array with the new column values.
        """
        
        if tablenum is None:
            for i in range(self.numtab):
                new_table = sd_fits_utils.update_table_column(self.table[i], 
                                                              column_name, 
                                                              column_vals)
                self.table[i] = new_table
        else:
            new_table = sd_fits_utils.update_table_column(self.table[tablenum], 
                                                          column_name, 
                                                          column_vals)
            self.table[tablenum] = new_table
        
            
