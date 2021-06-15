"""
Mapping functions.
"""

from groundhog.scan import Scan
from groundhog import sd_fits_utils


def map_with_ref(sdfits, scans, ref_scan, 
                 ifnum=0, plnum=0, fdnum=0, 
                 method='gbtidl', avgf_min=256):
    """
    Mapping with a reference position.
    
    Parameters
    ----------
    sdfits :
        
    scans : list
        Scan numbers.
    ref_scan : `Scan` object
        Reference scan.
    ifnum : int, optional
        Spectral window number.
    plnum : int, optional
        Polarization number.
    fdnum : int, optional
        Beam number.
    method : {'vector', 'gbtidl'}, optional
        Method used to compute the source temperature.
        If set to ``'vector'`` it will use Eq. (16) of
        Winkel et al. (2012). If set to ``'gbtidl'`` it
        will use the same method as GBTIDL.
        The default is ``'gbtidl'``.
    avgf_min : int, optional
        Averaging factor for the reference scan.
        The reference spectrum will be averaged by this
        factor in frequency.
    
    Returns
    -------
    cal_scans : `Scan` object
        Calibrated scans.
    """
    
    # Prepare variables from reference scan.
    ref = ref_scan.table['DATA']
    ref_tsys = ref_scan.table['TSYS']
    ref_exposure = ref_scan.table['EXPOSURE']
    
    # Prepare the mapping scans.
    sig_on = sdfits.get_scans(scans, ifnum=ifnum, plnum=plnum, fdnum=fdnum, cal='T')
    sig_off = sdfits.get_scans(scans, ifnum=ifnum, plnum=plnum, fdnum=fdnum, cal='F')
    
    # Add the spectra with the noise diode on and off.
    sig = (sig_on.data + sig_off.data)/2.
    sig_exposure = sig_on.table['EXPOSURE'] + sig_off.table['EXPOSURE']

    # Calibrate the mapping scans.
    ta = ref_tsys*(sig - ref)/ref
    exp = sig_exposure*ref_exposure / (sig_exposure + ref_exposure)

    # Create a new table with the calibrated data.
    cal_table = sig_on.table.copy()
    cal_table['EXPOSURE'] = exp
    cal_table['TSYS'] = ref_tsys
    cal_table['DATA'] = ta
    cal_scans = Scan(cal_table)
    
    return cal_scans
    
