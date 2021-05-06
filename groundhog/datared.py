"""
Data reduction functions.
"""

import warnings
import numpy as np

from groundhog import utils
from groundhog import spectral_axis
from groundhog.fluxscales import calibrators


def gbtidl_sigref2ta(sig, ref, tsys):
    """
    """
    
    return tsys[:,np.newaxis]*(sig - ref)/ref


def gbtidl_tsys(ref_on, ref_off, tcal):
    """
    """
    
    nchan = ref_on.shape[1]
    ch0 = int(nchan*0.1)
    chf = -int(nchan*0.1) + 1 # Python indexing is exclusive, IDL inclusive.
    
    return tcal*np.average(ref_off[:,ch0:chf], axis=1)/np.average(ref_on[:,ch0:chf] - ref_off[:,ch0:chf], axis=1) + tcal/2.
    

def get_ps(sdfits, scan, ifnum=0, intnum=None, plnum=0, method='freqdep'):
    """
    
    Parameters
    ----------
    sdfits :
        
    scan : int
        Scan number.
    plnum : int
        Polarization number.
    method : {'freqdep', 'classic'}, optional
        Method used to compute the source temperature.
        If set to ``'freqdep'`` it will use Eq. (16) of
        Winkel et al. (2012). If set to ``'classic'`` it
        will use the same method as GBTIDL.
        The default is ``'freqdep'``.
    
    Returns
    -------
    
    """

    ps_scan = sdfits.get_scans(scan, ifnum=ifnum, intnum=intnum, plnum=plnum)
    rows = ps_scan.table
    obsmode = rows["OBSMODE"]
    last_on = rows["LASTON"]
    last_off = rows["LASTOFF"]
    procnum = rows["PROCSEQN"]
    source = np.unique(rows['OBJECT'])[0]
    tcal = np.average(rows['TCAL'], axis=0)
    procname, swstate, swtchsig = obsmode[0].split(':')
    
    if procname not in ["OffOn", "OnOff"]:
        warnings.warn("Selected scan is not OnOff or OffOn, it is: {}"
                      "Cannot get Tcal from this scan.".format(procname))
        return None
    
    scan_on, scan_off = utils.get_ps_scan_pair(scan, procnum, procname)
    
    sou_on = sdfits.get_scans(scan_on, sig="T", cal="T", ifnum=ifnum, intnum=intnum, plnum=plnum)
    sou_off = sdfits.get_scans(scan_on, sig="T", cal="F", ifnum=ifnum, intnum=intnum, plnum=plnum)
    off_on = sdfits.get_scans(scan_off, sig="T", cal="T", ifnum=ifnum, intnum=intnum, plnum=plnum)
    off_off = sdfits.get_scans(scan_off, sig="T", cal="F", ifnum=ifnum, intnum=intnum, plnum=plnum)
    
    if method == 'freqdep':
        
        sou_on.average()
        sou_off.average()
        off_on.average()
        off_off.average()
    
        # Compute the kappa factor (Winkel et al. 2012).
        kappa_off = np.ma.power(off_on.data/off_off.data - 1., -1.)
        
        # Compute the source temperature (Eq. (16) in Winkel et al. 2012).
        tsou_on = (kappa_off + 1.)*tcal*(sou_on.data - off_on.data)/off_on.data
        tsou_off = kappa_off*tcal*(sou_off.data - off_off.data)/off_off.data
        # Average.
        tsou = 0.5*(tsou_on + tsou_off)
        
    elif method == 'classic':
        
        # Eqs. (1) and (2) from Braatz (2009, GBTIDL calibration guide). 
        tsys = gbtidl_tsys(off_on.data, off_off.data, tcal)
        sig = 0.5*(sou_on.data + sou_off.data)
        ref = 0.5*(off_on.data + off_off.data)
        tsou_int = gbtidl_sigref2ta(sig, ref, tsys)
        tsig_sou = 0.5*(sou_on.table["EXPOSURE"] + sou_off.table["EXPOSURE"])
        tsig_off = 0.5*(off_on.table["EXPOSURE"] + off_off.table["EXPOSURE"])
        tsig = 0.5*(tsig_sou + tsig_off)
        dnu = np.mean(sou_on.table["CDELT1"])
        tsou = np.average(tsou_int, axis=0, weights=dnu*tsig*np.power(tsys, -2.))
    
    return tsou


def get_tcal(sdfits, scan, ifnum=0, intnum=None, plnum=0, scale="Perley-Butler 2017", units="K"):
    """
    """
    
    cal_scan = sdfits.get_scans(scan, ifnum=ifnum, intnum=intnum, plnum=plnum)
    cal_rows = cal_scan.table
    obsmode = cal_rows["OBSMODE"]
    last_on = cal_rows["LASTON"]
    last_off = cal_rows["LASTOFF"]
    procnum = cal_rows["PROCSEQN"]
    source = np.unique(cal_rows['OBJECT'])[0]
    
    procname, swstate, swtchsig = obsmode[0].split(':')
    
    if procname not in ["OffOn", "OnOff"]:
        warnings.warn("Selected scan is not OnOff or OffOn, it is: {}"
                      "Cannot get Tcal from this scan.".format(procname))
        return None
    
    if "PSWITCH" not in swstate:
        warnings.warn("Selected scan is not position switched "
                      "check results.")
    
    scan_on, scan_off = utils.get_ps_scan_pair(scan, procnum, procname)
    
    # Get the On and Off source scans with the noise diode On and Off.
    sou_on = sdfits.get_scans(scan_on, sig="T", cal="T", ifnum=ifnum, intnum=intnum, plnum=plnum)
    sou_off = sdfits.get_scans(scan_on, sig="T", cal="F", ifnum=ifnum, intnum=intnum, plnum=plnum)
    off_on = sdfits.get_scans(scan_off, sig="T", cal="T", ifnum=ifnum, intnum=intnum, plnum=plnum)
    off_off = sdfits.get_scans(scan_off, sig="T", cal="F", ifnum=ifnum, intnum=intnum, plnum=plnum)
    
    sou_on.average()
    sou_off.average()
    off_on.average()
    off_off.average()
    
    # Compute the kappa factor (Winkel et al. 2012).
    kappa_off = np.ma.power(off_on.data/off_off.data - 1., -1.)
    
    ta_sou_on = calibrators.compute_sed(sou_on.freq, scale, source, units=units)
    ta_sou_off = calibrators.compute_sed(off_on.freq, scale, source, units=units)
    
    # Compute the temperature of the noise diode (Eq. (76) in Winkel et al. 2012).
    # Using the observations with the noise diode off.
    tcal_off = ta_sou_off/(kappa_off*(sou_off.data - off_off.data)/off_off.data)
    # Using the observations with the noise diode on.
    tcal_on = ta_sou_on/((kappa_off + 1.)*(sou_on.data - off_on.data)/off_on.data)
    # Average the results.
    tcal = (tcal_off/np.ma.std(tcal_off)**2. + tcal_on/np.ma.std(tcal_on)**2.) / \
           (1./np.ma.std(tcal_off)**2. + 1./np.ma.std(tcal_on)**2.)

    return tcal
        
        
    
    
    
