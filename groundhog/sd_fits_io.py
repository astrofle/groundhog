
"""
Single Dish FITS I/O methods.
https://fits.gsfc.nasa.gov/registry/sdfits.html
"""

from astropy.io import fits

from groundhog.sd_fits import SDFITS


def read_sdfits(filename, ext='SINGLE DISH'):
    """
    """
    
    #hdu = fits.open(filename)
    #table = hdu[ext].data
    #head = hdu[ext].header
    table, head = _read_sdfits(filename, ext=ext)
    sdfits = SDFITS(table, head)
    
    return sdfits


def _read_sdfits(filename, ext='SINGLE DISH'):
    """
    """
    
    hdu = fits.open(filename)
    table = hdu[ext].data
    head = hdu[ext].header
    
    return table, head
    

def write_sdfits(filename, sdfits, overwrite=False):
    """
    """
    
    table = sdfits.table
    header = sdfits.header
    
    _write_sdfits(filename, table, header, overwrite=overwrite)
    
    
def _write_sdfits(filename, table, header, overwrite=False):
    """
    """
    
    new_hdu = fits.BinTableHDU.from_columns(table)
    
    for k,v in header.items(): 
        if k not in new_hdu.header.keys(): 
            new_hdu.header[k] = v
            
    new_hdu.writeto(filename, overwrite=overwrite)
