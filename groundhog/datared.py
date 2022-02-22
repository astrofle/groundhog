"""
Data reduction functions.
"""

import warnings
import numpy as np

from groundhog import utils
from groundhog.scan import Scan
from groundhog import sd_fits_utils
from groundhog import spectral_axis
from groundhog.fluxscales import calibrators


def gbtidl_sigref2ta(sig, ref, tsys):
    """
    """
    
    return tsys[:,np.newaxis]*(sig - ref)/ref


def gbtidl_tsys(ref_on, ref_off, tcal, use80=True):
    """
    """
    
    nchan = ref_on.shape[1]
    ch0 = int(nchan*0.1)
    chf = -int(nchan*0.1) + 1 # Python indexing is exclusive, IDL inclusive.
    rng = slice(ch0,chf,1)

    if not use80:
        rng = slice(0,None,1)

    return tcal*np.ma.average(ref_off[:,rng], axis=1)/np.average(ref_on[:,rng] - ref_off[:,rng], axis=1) + tcal/2.


def classic_tsys(ref_on, ref_off, tcal):
    """
    """
    
    nchan = ref_on.shape[1]
    ch0 = int(nchan*0.1)
    chf = -int(nchan*0.1) + 1 # Python indexing is exclusive.
    
    tsys_tcal = (ref_on + ref_off - np.average((ref_on + ref_off)[:,ch0:chf], axis=1)[:,np.newaxis])/(2.*np.average((ref_on + ref_off)[:,ch0:chf], axis=1)[:,np.newaxis])
    
    return tcal*np.average(tsys_tcal, axis=1)


def get_kappa(tcal_on, tcal_off, avgf=1):
    """
    Computes the kappa factor (Eq. (14) in Winkel et al. 2012).
    
    Parameters
    ----------
    tcal_on : array
        Power with the noise diode on.
    tcal_off : array
        Power with the noise diode off.
    avgf : int
        Average the ratio `tcal_on/tcal_off` by this amount.
    
    Returns
    -------
    kappa : array
        Kappa factor as defined by Eq. (14) in Winkel et al. 2012.
    """
    
    nchan = len(tcal_on)
    
    # Compute the kappa factor (Winkel et al. 2012).
    off_ratio = tcal_on/tcal_off
    # Average in frequency to increase the SNR.
    off_ratio = off_ratio.reshape(nchan//avgf, avgf).mean(axis=1)
    kappa = np.ma.power(off_ratio - 1., -1.)  
    
    return kappa
    

def get_ps(sdfits, scan, ifnum=0, intnum=None, plnum=0, fdnum=0, method='vector', avgf_min=256):
    """
    
    Parameters
    ----------
    sdfits :
        
    scan : int
        Scan number.
    plnum : int
        Polarization number.
    method : {'vector', 'classic'}, optional
        Method used to compute the source temperature.
        If set to ``'vector'`` it will use Eq. (16) of
        Winkel et al. (2012). If set to ``'classic'`` it
        will use the same method as GBTIDL.
        The default is ``'vector'``.
    
    Returns
    -------
    
    """

    ps_scan = sdfits.get_scans(scan, ifnum=ifnum, intnum=intnum, plnum=plnum)
    rows = ps_scan.array
    obsmode = rows["OBSMODE"]
    last_on = rows["LASTON"]
    last_off = rows["LASTOFF"]
    procnum = rows["PROCSEQN"]
    source = np.unique(rows['OBJECT'])[0]
    tcal = np.average(rows['TCAL'], axis=0)
    procname, swstate, swtchsig = obsmode[0].split(':')
    
    if procname not in ["OffOn", "OnOff"]:
        warnings.warn(f"Selected scan is not OnOff or OffOn, it is: {procname}."
                      f"Cannot get Tcal from this scan.")
        return None
    
    scan_on, scan_off = utils.get_ps_scan_pair(scan, procnum, procname)
    
    sou_on = sdfits.get_scans(scan_on, sig="T", cal="T", ifnum=ifnum, intnum=intnum, plnum=plnum)
    sou_off = sdfits.get_scans(scan_on, sig="T", cal="F", ifnum=ifnum, intnum=intnum, plnum=plnum)
    off_on = sdfits.get_scans(scan_off, sig="T", cal="T", ifnum=ifnum, intnum=intnum, plnum=plnum)
    off_off = sdfits.get_scans(scan_off, sig="T", cal="F", ifnum=ifnum, intnum=intnum, plnum=plnum)
    
    if method == 'vector':
        
        sou_on.average()
        sou_off.average()
        off_on.average()
        off_off.average()
        off_freq = off_off.freq
        sou_freq = sou_on.freq
        
        nchan = off_on.data.shape[0]
        facs = utils.factors(nchan)
        avgf = np.min(facs[facs >= avgf_min])
        
        kappa_off = get_kappa(off_on.data, off_off.data, avgf=avgf)
        kappa_freq = off_freq.reshape(nchan//avgf, avgf).mean(axis=1)
        
        # Interpolate back to high frequency resolution.
        pt = np.argsort(kappa_freq)
        pi = np.argsort(sou_freq)
        kappa_interp = np.interp(sou_freq.to('Hz').value[pi], kappa_freq.to('Hz').value[pt], kappa_off)
        
        # Compute the source temperature (Eq. (16) in Winkel et al. 2012).
        tsou_on = (kappa_interp + 1.)*tcal*(sou_on.data - off_on.data)/off_on.data
        tsou_off = kappa_interp*tcal*(sou_off.data - off_off.data)/off_off.data
        # Average.
        tsou = 0.5*(tsou_on + tsou_off)
        
    elif method == 'gbtidl':
        
        # Eqs. (1) and (2) from Braatz (2009, GBTIDL calibration guide)
        # https://www.gb.nrao.edu/GBT/DA/gbtidl/gbtidl_calibration.pdf
        tsys = gbtidl_tsys(off_on.data, off_off.data, tcal)
        sig = 0.5*(sou_on.data + sou_off.data)
        ref = 0.5*(off_on.data + off_off.data)
        ta = gbtidl_sigref2ta(sig, ref, tsys)
        tint_sou = 0.5*(sou_on.array["EXPOSURE"] + sou_off.array["EXPOSURE"])
        tint_off = 0.5*(off_on.array["EXPOSURE"] + off_off.array["EXPOSURE"])
        tint = 0.5*(tint_sou + tint_off)
        dnu = np.mean(sou_on.array["CDELT1"])
        tsou = np.average(ta, axis=0, weights=dnu*tint*np.power(tsys, -2.))
        
    elif method == 'classic':
        tsys = classic_tsys(off_on.data, off_off.data, tcal)
        ta_on = (sou_on.data - off_on.data)/off_on.data*(tsys[:,np.newaxis] + tcal)
        ta_off = (sou_off.data - off_off.data)/off_off.data*(tsys[:,np.newaxis])
        tint_sou = 0.5*(sou_on.array["EXPOSURE"] + sou_off.array["EXPOSURE"])
        tint_off = 0.5*(off_on.array["EXPOSURE"] + off_off.array["EXPOSURE"])
        tint = 0.5*(tint_sou + tint_off)
        dnu = np.mean(sou_on.array["CDELT1"])
        ta_on = np.average(ta_on, axis=0, weights=dnu*tint_sou*np.power(tsys, -2.))
        ta_off = np.average(ta_off, axis=0, weights=dnu*tint_off*np.power(tsys, -2.))
        tsou = 0.5*(ta_on + ta_off)
        
    return tsou


def get_tcal(sdfits, scan, ifnum=0, intnum=None, plnum=0, scale="Perley-Butler 2017", units="K",
             avgf_min=16):
    """
    """
    
    cal_scan = sdfits.get_scans(scan, ifnum=ifnum, intnum=intnum, plnum=plnum)
    cal_rows = cal_scan.array
    obsmode = cal_rows["OBSMODE"]
    last_on = cal_rows["LASTON"]
    last_off = cal_rows["LASTOFF"]
    procnum = cal_rows["PROCSEQN"]
    source = np.unique(cal_rows['OBJECT'])[0].strip()
    
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
    off_off.get_freq()
    sou_on.get_freq()
    off_freq = off_off.freq[0]
    sou_freq = sou_on.freq[0]
    
    nchan = off_on.data.shape[0]
    facs = utils.factors(nchan)
    avgf = np.min(facs[facs >= avgf_min])
    
    kappa_off = get_kappa(off_on.data, off_off.data, avgf=avgf)
    kappa_freq = off_freq.reshape(nchan//avgf, avgf).mean(axis=1)
    
    # Interpolate back to high frequency resolution.
    pt = np.argsort(kappa_freq)
    pi = np.argsort(sou_freq)
    kappa_interp = np.interp(sou_freq.to('Hz').value[pi], kappa_freq.to('Hz').value[pt], kappa_off)
    
    ta_sou_on = calibrators.compute_sed(sou_freq, scale, source, units=units)
    ta_sou_off = calibrators.compute_sed(off_freq, scale, source, units=units)
    
    # Compute the temperature of the noise diode (Eq. (76) in Winkel et al. 2012).
    # Using the observations with the noise diode off.
    tcal_off = ta_sou_off/(kappa_interp*(sou_off.data - off_off.data)/off_off.data)
    # Using the observations with the noise diode on.
    tcal_on = ta_sou_on/((kappa_interp + 1.)*(sou_on.data - off_on.data)/off_on.data)
    # Average the results.
    tcal = (tcal_off/np.ma.std(tcal_off)**2. + tcal_on/np.ma.std(tcal_on)**2.) / \
           (1./np.ma.std(tcal_off)**2. + 1./np.ma.std(tcal_on)**2.)

    return tcal
        
        
def prepare_mapping_off(sdfits, scan, ifnum=0, intnum=None, plnum=0, fdnum=0):
    """
    """
    
    ref_on = sdfits.get_scans(scan, ifnum=ifnum, plnum=plnum, fdnum=fdnum, cal='T')
    ref_off = sdfits.get_scans(scan, ifnum=ifnum, plnum=plnum, fdnum=fdnum, cal='F')

    tcal = np.average(ref_on.table['TCAL'], axis=0)

    ref_add = (ref_on.data + ref_off.data)/2.
    ref_exposures = ref_on.table['EXPOSURE'] + ref_off.table['EXPOSURE']

    tsys_ref = gbtidl_tsys(ref_on.data, ref_off.data, tcal)

    weight = ref_exposures / tsys_ref**2
    tsys = np.average(tsys_ref, weights=weight)
    ref = np.ma.average(ref_add, weights=weight, axis=0)
    ref_exposure = ref_exposures.sum()
    
    # Create a scan object for the reference position
    # using the average values from the integrations.
    ref_table = ref_on.table.copy()[0]
    ref_table.setfield('DATA', ref)
    ref_table.setfield('TSYS', tsys)
    ref_table.setfield('EXPOSURE', ref_exposure)
    for k in np.arange(1,4,1):
        avg_col = np.array([np.average(ref_on.table.field(f'crval{k}'), weights=weight)])
        ref_table.setfield(f'crval{k}', avg_col)
    ref_scan = Scan(ref_table)
    ref_scan.get_freq()
    
    return ref_scan
    
