
"""
Utility functions to handle SDFITS files.
"""


import numpy as np
import numpy.lib.recfunctions as rfn

from collections import namedtuple

from astropy.io import fits


def build_index(hdu, extname='SINGLE_DISH', max_chunk=10000):
    """
    """

    keys = ['scan', 'cal', 'plnum', 'sig', 'ifnum', 'fdnum',
            'object', 'velocity', 'obsmode', 'procseqn', 'restfreq',
            'azimuth', 'elevatio']
    dtypes = [int, str, int, str, int, int, 
              "<U10", float, "<U20", int, float, 
              float, float]
    
    index = {}

    for i,h in enumerate(hdu):
        
        if h.get_extname() == extname:
            
            index[i] = {}
            nrows = h.get_info()['nrows']

            # Create empty arrays to store the values.
            for j,k in enumerate(keys):
                index[i][k] = np.empty(nrows, dtype=dtypes[j])
          
            # Extract the contents of the HDU in chunks.
            for c in range(0, int(np.ceil(nrows/max_chunk))+1):
                c0 = c*max_chunk
                cf = (c+1)*max_chunk

                # Get the values from the HDU.
                cols = h[keys][c0:cf]
                for k in keys:
                    #cols = h[k][c0:cf]
                    index[i][k][c0:cf] = cols[k.upper()]
                    #index[i][k][c0:cf] = h[k][c0:cf]

    # Set unique values from the rows,
    # remove extra characters from the object name,
    # and set the rest frequency to GHz.
    for i,h in enumerate(hdu):
        if h.get_extname() == extname:
            for j,k in enumerate(keys):
                index[i][f"u{k}"] = np.array(list(set(index[i][k])), dtype=dtypes[j])

            index[i]['object'] = np.array([o.strip() for o in index[i]['object']])
            index[i]['restfreq'] *= 1e-9 # GHz

    return index


def get_sdfits_ext(sdfits, extname="SINGLE_DISH"):
    """
    """

    if len(sdfits.hdu) == 2:
        extnum = 1
    elif len(sdfits.hdu) > 2:
        for i in range(1,len(sdfits.hdu)):
            if sdfits.hdu[i].get_extname() == extname:
                extnum = i
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

    return extnum


def get_rows(index, scans=None, ifnum=None, sig=None, cal=None, plnum=None, fdnum=None):
    """
    """

    rows = np.arange(len(index['scan']))

    if scans is None and ifnum is None and sig is None and cal is None and plnum is None:
        mask = False

    else:
        if scans is not None:
            mask = np.isin(index["scan"], scans)
        else:
            scans = index["uscan"]
            mask = np.isin(index["scan"], scans)
        
        if ifnum is not None:
            mask = mask & (index["ifnum"] == ifnum)
            
        if sig is not None:
            mask = mask & (index["sig"] == sig)
            
        if cal is not None:
            mask = mask & (index["cal"] == cal)

        if plnum is not None:
            mask = mask & np.isin(index["plnum"], plnum)

        if fdnum is not None:
            mask = mask & np.isin(index["fdnum"], fdnum)

    return rows[mask]



def get_table_mask(table, scans=None, ifnum=None, sig=None, cal=None, plnum=None, fdnum=None):
    """
    """
    
    if scans is None and ifnum is None and sig is None and cal is None and plnum is None:
        mask = None
        
    else:
        if scans is not None:
            mask = np.isin(table["SCAN"], scans)
        else:
            scans = np.unique(table["SCAN"])
            mask = np.isin(table["SCAN"], scans)
        
        if ifnum is not None:
            mask = mask & (table["IFNUM"] == ifnum)
            
        if sig is not None:
            mask = mask & (table["SIG"] == sig)
            
        if cal is not None:
            mask = mask & (table["CAL"] == cal)
            
        if plnum is not None:
            mask = mask & np.isin(table["PLNUM"], plnum)
            
        if fdnum is not None:
            mask = mask & np.isin(table["FDNUM"], fdnum)
    
    return mask


def update_column_fmt(column, new_array):
    """
    """

    fmt = column.format

    if fmt[-1] == 'A':
        new_fmt = fmt
    elif len(new_array.shape) == 1:
        new_fmt = '1{}'.format(fmt[-1])
    else:
        new_fmt = '{}{}'.format(new_array.shape[1], fmt[-1])

    return new_fmt


def update_table_column(table, column, new_array):
    """
    Parameters
    ----------
    table : `astropy.io.fits.fitsrec`
        Contents of the SDFITS file to update.
    column : str
        Name of the column to update.
    new_array : np.ndarray
        Array with the new column values.
    """
    
    # Get the original table values.
    cols = table.columns
   
    # Define the fmt of the updated column.
    fmt = update_column_fmt(cols[column], new_array)
 
    # Delete the column we will update.
    cols.del_col(column)
    
    # Define the updated column and add it to the other columns.
    new_col = fits.Column(name=column, format=fmt, array=new_array)
    cols.add_col(new_col)
    
    # Make it a new table.
    new_hdu = fits.BinTableHDU.from_columns(cols)
    new_table = new_hdu.data
    
    return new_table


def append_table_column(table, column, data, dtype, loc=None):
    """
    """

    cols_dt_list = table.dtype.descr

    # Do not try to update existing columns.
    if column in cols_dt_list:
        print("Column already present.")
        print("Will not overwrite.")
        return table

    # Define the data type for the new table with the 
    # new column in the requested location.
    if loc is None:
        loc = len(cols_dt_list)
    cols_dt_list.insert(loc, dtype)
    new_dt = np.dtype(cols_dt_list)

    # Create a new table with the same number of rows.
    new_table = np.empty(len(table), dtype=new_dt)
    
    # Copy the contents from the original table into the new table.
    for n in table.dtype.names:
        new_table[n] = table[n]
    new_table[column] = data

    return new_table


def split_channel_range(table, ch0, chf, dch=1):
    """
    """

    # Find the number of channels and 
    # define how many channels the selection will have.
    nrow, nchan = table['DATA'].shape
    fslice = slice(ch0, chf, dch)
    chan_slice = fslice.indices(nchan)
    nchan_sel = (chan_slice[1] - chan_slice[0])//chan_slice[2]

    # Remove the DATA column from the table.
    nodata_table = rfn.drop_fields(table, 'DATA')
    # Copy the column definitions as a list.
    nodata_table_dt = nodata_table.dtype.descr
    # Concatenate the column definitions with the new data shape.
    new_dt = np.dtype(nodata_table_dt[:6] + [('DATA', '>f4', (nchan_sel,))] + nodata_table_dt[6:])
    # Create a new table with the same number of rows.
    new_table = np.empty(nrow, dtype=new_dt)

    # Fill the new table with the old contents, 
    # and the DATA selection.
    for n in nodata_table.dtype.names:
        new_table[n] = nodata_table[n]
    new_table['DATA'] = table['DATA'][:,fslice]
    # Update the frequency axis and bandwidth.
    new_table['CRPIX1'] -= ch0
    new_table['BANDWID'] = new_table['FREQRES'] * nchan_sel

    return new_table

