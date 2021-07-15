
"""
Utility functions to handle SDFITS files.
"""


import numpy as np

from collections import namedtuple

from astropy.io import fits


def build_index(hdu, ext='SINGLE DISH'):
    """
    """
    
    keys = ['scan', 'cal', 'plnum', 'sig', 'ifnum', 'fdnum']
    
    index = {}

    for i,h in enumerate(hdu):
        if h.name == ext:
            index[i] = {}
            for k in keys:
                index[i][k] = h.data.field(k)
                index[i][f"u{k}"] = np.array(list(set(h.data.field(k))))
    
    return index


def parse_sdfits(table):
    """
    
    """
    
    uscans = np.array(list(set(table['SCAN'])), dtype=int)
    ucal = np.array(list(set(table['CAL'])), dtype=str)
    uplnum = np.array(list(set(table['PLNUM'])), dtype=int)
    usig = np.array(list(set(table['SIG'])), dtype=str)
    uifnum = np.array(list(set(table['IFNUM'])), dtype=str)
    ufdnum = np.array(list(set(table['FDNUM'])), dtype=str)
    
    
    Unique = namedtuple('UniqueEntries', ['scan', 'cal', 'plnum', 'sig', 'ifnum', 'fdnum'])
    uniques = Unique(scan=uscans, cal=ucal, plnum=uplnum, sig=usig, ifnum=uifnum, fdnum=ufdnum)
    
    return uniques


def make_summary(table):
    """
    Creates a summary of the contents of an SDFITS table.
    It tries to imitate the `summary` function in GBTIDL.
    
    Parameters
    ----------
    table : `astropy.io.fits.fitsrec`
        Contents of the SDFITS file to summarize.
        
    Returns
    -------
    summary : `Summary` object
        Summary(Scan, Source, Vel, Proc, Seq, RestF, nIF, nInt, nFd, Az, El)
    """
    
    uniques = parse_sdfits(table)
    
    Summary = namedtuple('Summary', ['Scan', 'Source', 'Vel', 'Proc', 'Seq',
                                     'RestF', 'nIF', 'nInt', 'nFd', 'Az', 'El'])
    summary = Summary(Scan=[], Source=[], Vel=[], Proc=[], Seq=[],
                      RestF=[], nIF=[], nInt=[], nFd=[], Az=[], El=[])
    
    for scan in uniques.scan:
        
        mask = (table['SCAN'] == scan)
        summary.Scan.append(scan)
        summary.Source.append(np.unique(table['OBJECT'][mask])[0])
        summary.Vel.append(np.mean(table['VELOCITY'][mask]))
        summary.Proc.append(np.unique(table['OBSMODE'][mask]).split(':')[0][0])
        summary.Seq.append(np.unique(table['PROCSEQN'][mask])[0])
        summary.RestF.append(np.unique(table['RESTFREQ'][mask])[0])
        summary.nIF.append(len(np.unique(table['IFNUM'][mask])))
        summary.nInt.append(mask.sum())
        summary.nFd.append(np.unique(table['FDNUM'][mask])[0])
        summary.Az.append(np.mean(table['AZIMUTH'][mask]))
        summary.El.append(np.mean(table['ELEVATIO'][mask]))
        
    return summary


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
    fmt = cols[column].format
    
    # Delete the column we will update.
    cols.del_col(column)
    
    if fmt[-1] == 'A':
        new_fmt = fmt
    elif len(new_array.shape) == 1:
        new_fmt = '1{}'.format(fmt[-1])
    else:
        new_fmt = '{}{}'.format(new_array.shape[1], fmt[-1])
    
    # Define the updated column.
    new_col = fits.Column(name=column, format=new_fmt, array=new_array)
    
    cols.add_col(new_col)
    
    # Make it a new table.
    new_hdu = fits.BinTableHDU.from_columns(cols)
    new_table = new_hdu.data
    
    return new_table
