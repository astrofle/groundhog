
"""
"""

import numpy as np
import numpy.lib.recfunctions as rfn

from collections import namedtuple

import fitsio
import pandas as pd
from astropy.io import fits

from groundhog.scan import Scan
from groundhog import datared
from groundhog import sd_fits_utils


class SDFITS:
    
    __name__ = "SDFITS"
    
    def __init__(self, hdu=None, index=None, header=None):
        
        self.hdu = hdu
        self.index = index
        self._summary = None
       
    
    def get_rows(self, extnum, scans=None, ifnum=None, sig=None, cal=None, plnum=None, fdnum=None, intnum=None):
        """
        """

        rows = sd_fits_utils.get_rows(self.index[extnum],
                                      scans=scans, ifnum=ifnum,
                                      sig=sig, cal=cal,
                                      plnum=plnum, fdnum=fdnum)

        table = self.hdu[extnum][rows]
        #index = self.index[extnum][rows]
        
        #outfits = fitsio.FITS()

        return table


    def load(self, filename, max_chunk=10000, mode='r', extname='SINGLE_DISH'):
        """
        """
        
        self.hdu = fitsio.FITS(filename, mode)
        self.index = sd_fits_utils.build_index(self.hdu, extname=extname, 
                                               max_chunk=max_chunk)

    
    def make_summary(self):
        """
        Creates a summary of the contents of an SDFITS table.
        It tries to imitate the `summary` function in GBTIDL.
        """

        Summary = namedtuple('Summary', ['Scan', 'Source', 'Vel', 'Proc', 'Seq',
                                         'RestF', 'nIF', 'nInt', 'nFd', 'Az', 'El'])
        summary = Summary(Scan=[], Source=[], Vel=[], Proc=[], Seq=[],
                          RestF=[], nIF=[], nInt=[], nFd=[], Az=[], El=[])

        for k in self.index.keys():
            index = self.index[k]
            for i,scan in enumerate(index['uscan']):
                summary.Scan.append(scan)
                mask = (index['scan'] == scan)
                summary.Source.append(list(set(index['object'][mask]))[0])
                summary.Vel.append(np.mean(index['velocity'][mask]))
                summary.Proc.append(list(set(index['obsmode'][mask]))[0].split(':')[0])
                summary.Seq.append(list(set(index['procseqn'][mask]))[0])
                summary.RestF.append(list(set(index['restfreq'][mask]))[0])
                summary.nIF.append(len(list(set(index['ifnum'][mask]))))
                # Need to fix the number of integrations. Currently it reports #rows/#spws.
                summary.nInt.append(mask.sum()/summary.nIF[i])
                summary.nFd.append(list(set(index['fdnum'][mask]))[0])
                summary.Az.append(np.mean(index['azimuth'][mask]))
                summary.El.append(np.mean(index['elevatio'][mask]))

        self._summary = summary


    def get_channels(self, ch0, chf, dch=1, extname='SINGLE_DISH'):
        """
        Returns the contents of the SDFITS file with its DATA containing only the
        selected channels.
        Only dch=1 is supported at the moment.
        This will onlyreturn the contents of the first extension matching extname.
        """

        if dch > 1 or dch <= 0:
            print("Not implemented.")
            return

        # Find the extension number and read its contents. 
        extnum = sd_fits_utils.get_sdfits_ext(self, extname=extname)
        table = self.hdu[extnum][:]
        head = self.hdu[extnum].read_header()

        # Find the number of channels and 
        # define how many channels the selection will have.
        #nchan = int(head['TFORM7'][:-1])
        nrow, nchan = table['DATA'].shape
        fslice = slice(ch0, chf, dch)
        chan_slice = fslice.indices(nchan)
        nchan_sel = (chan_slice[1] - chan_slice[0])//chan_slice[2]

        # Remove the DATA column from the table.
        nodata_cols = rfn.drop_fields(table, 'DATA')
        # Copy the column definitions as a list.
        nodata_cols_dt = nodata_cols.dtype.descr
        # Concatenate the column definitions with the new data shape.
        new_dt = np.dtype(nodata_cols_dt[:6] + [('DATA', '>f4', (nchan_sel,))] + nodata_cols_dt[6:])
        # Create a new table with the same number of rows.
        new_table = np.empty(head['NAXIS2'], dtype=new_dt)
        
        # Fill the new table with the old contents, 
        # and the DATA selection.
        for n in nodata_cols.dtype.names:
            new_table[n] = nodata_cols[n]
        new_table['DATA'] = table['DATA'][:,fslice]
        # Update the frequency axis and bandwidth.
        new_table['CRPIX1'] -= ch0
        new_table['BANDWID'] = new_table['FREQRES'] * nchan_sel

        return new_table


    def get_scans(self, scans, ifnum=None, sig=None, cal=None, plnum=None, fdnum=None, intnum=None):
        """
        Returns the rows in the SDFITS table for the requested scan numbers.
        
        Parameters
        ----------
        scans : int or array_like
            Scans to retrieve.
        
        Returns
        -------
        scan : `groundhog.Scan`
            Rows with the selected scans.
        """
        
        table = None
        extnum = None
        rows = None
        
        # Find the extension.
        # There's always the primary HDU, plus the table HDUs.
        #get_sdfits_ext(sdfits
        if len(self.hdu) == 2:
            extnum = 1
        elif len(self.hdu) > 2:
            for i in range(1,len(self.hdu)):
                if np.isin(self.index[i]['uscan'], scans).sum() > 0:
                    # If something changed during an observation, e.g.,
                    # number of channels, the data is put in separate
                    # tables in the same sdfits file. In principle,
                    # we could work with data that spans multiple
                    # configurations, but it is better to leave that to 
                    # the users to avoid making decisions for them.
                    if extnum is not None:
                        print("Scans span multiple configurations.")
                        print("This is not supported.")
                        return
                    extnum = i
        
        rows = sd_fits_utils.get_rows(self.index[extnum], 
                                      scans=scans, ifnum=ifnum, 
                                      sig=sig, cal=cal, 
                                      plnum=plnum, fdnum=fdnum)
        
        out_scans = self.hdu[extnum][rows]
        
        if intnum is not None:
            if not hasattr(intnum, "__len__"):
                intnum = [intnum]
            out_scans = out_scans[intnum]
        
        scan = Scan(out_scans)
        
        return scan
    
    
    def remove_edge_chans(self, hduidx=1, frac=0.2, chan0=None, chanf=None, dch=1):
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
        
        #for i,table in enumerate(self.hdu[1:]):
        table = self.hdu[hduidx]
        data = table["DATA"][:]
        nrow, nchan = data.shape #table['DATA'].shape

        #if len(data.shape) == 1:
        #    nchan = data.shape[0]
        #elif len(data.shape) == 2:
        #    nchan = data.shape[1]
            
        if chan0 is None:
            chan0 = int(nchan*frac/2)
        if chanf is None:
            chanf = int(nchan - nchan*frac/2)
        
        xslice = slice(chan0, chanf, dch)
        chan_slice = xslice.indices(nchan)
        nchan_sel = (chan_slice[1] - chan_slice[0])//chan_slice[2]

        # Remove the DATA column from the table.
        nodata_table = rfn.drop_fields(table[:], "DATA")
        # Copy the column definitions as a list.
        nodata_table_dt = nodata_table.dtype.descr
        # Concatenate the column definitions with the new data shape.
        new_dt = np.dtype(nodata_table_dt[:6] + [("DATA", '>f4', (nchan_sel,))] + nodata_table_dt[6:])
        # Create a new table with the same number of rows.
        new_table = np.empty(nrow, dtype=new_dt)

        # Fill the new table with the old contents, 
        # and the DATA selection.
        for n in nodata_table.dtype.names:
            new_table[n] = nodata_table[n]
        new_table["DATA"] = table["DATA"][:][:,xslice]
        # Update the frequency axis and bandwidth.
        new_table["CRPIX1"] -= chan0
        new_table["BANDWID"] = new_table["FREQRES"] * nchan_sel

        self.hdu[hduidx][:] = new_table

        ## Select only inner channels.
        #if len(data.shape) == 1:
        #    data = data[chan0:chanf]
        #elif len(data.shape) == 2:
        #    data = data[:,chan0:chanf]
        #
        ## Update DATA column.
        #new_table = sd_fits_utils.update_table_column(table, 'DATA', data)
        #self.hdu[i]['DATA'] = new_table
        #
        ## Update frequency axis header.
        #crp1 = table.field('crpix1')
        #crp1 -= chan0
        #new_table = sd_fits_utils.update_table_column(table, 'CRPIX1', crp1)
        #self.hdu[i]['DATA'] = new_table
            
    
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
        

    def summary(self):
        """
        """

        if self._summary is None:
            self.make_summary()

        pd.set_option('display.max_rows', None)
        pd.options.display.float_format = '{:,.3f}'.format
        pd.set_option('display.expand_frame_repr', False)
        df = pd.DataFrame(data=self._summary._asdict())

        print(df.to_string(index=False))
