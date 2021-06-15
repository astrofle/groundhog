
"""
Single Dish FITS I/O methods.
https://fits.gsfc.nasa.gov/registry/sdfits.html
"""

import numpy as np
from datetime import datetime

from astropy.io import fits
from astropy.time import Time

from groundhog.sd_fits import SDFITS


def make_sdfits_primary_header(date=None, telescope='NRAO_GBT', backend='VEGAS'):
    """
    """
    
    primary_sdfits_header =  "SIMPLE  =                    T / file does conform to FITS standard             BITPIX  =                    8 / number of bits per data pixel                  NAXIS   =                    0 / number of data axes                            EXTEND  =                    T / FITS dataset may contain extensions            COMMENT   FITS (Flexible Image Transport System) format is defined in 'AstronomyCOMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H DATE    = '2017-05-17T04:27:03' / date and time this HDU was created            ORIGIN  = 'NRAO Green Bank'    / origin of observation                          TELESCOP= 'NRAO_GBT'           / the telescope used                             INSTRUME= 'VEGAS   '           / backend                                        SDFITVER= 'sdfits ver1.21'     / this file was created by sdfits                FITSVER = '1.9     '           / FITS definition version                        END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "
    
    header = fits.Header.fromstring(primary_sdfits_header)
    
    if date is None:
        ut = Time(datetime.utcnow(), scale='utc')
        ut.format = 'fits'
        header['DATE'] = ut.value
    header['TELESCOP'] = telescope
    header['INSTRUME'] = backend
    
    return header


def make_sdfits_table_header():
    """
    """
    
    table_header = 'XTENSION= \'BINTABLE\'           / binary table extension                         BITPIX  =                    8 / 8-bit bytes                                    NAXIS   =                    2 / 2-dimensional binary table                     NAXIS1  =               131698 / width of table in bytes                        NAXIS2  =                15680 / number of rows in table                        PCOUNT  =                    0 / size of special data area                      GCOUNT  =                    1 / one data group (required keyword)              TFIELDS =                   74 / number of fields in each row                   COMMENT Start of SDFITS CORE keywords/columns.                                  TTYPE1  = \'OBJECT  \'           / name of source observed                        TFORM1  = \'32A     \'                                                            TUNIT1  = \'        \'                                                            TELESCOP= \'NRAO_GBT\'           / the telescope used                             TTYPE2  = \'BANDWID \'           / bandwidth                                      TFORM2  = \'1D      \'                                                            TUNIT2  = \'Hz      \'                                                            TTYPE3  = \'DATE-OBS\'           / date and time of observation start             TFORM3  = \'22A     \'                                                            TUNIT3  = \'        \'                                                            TTYPE4  = \'DURATION\'           / total integration duration in seconds          TFORM4  = \'1D      \'                                                            TUNIT4  = \'s       \'                                                            TTYPE5  = \'EXPOSURE\'           / effective int time (excludes blanking) in secs TFORM5  = \'1D      \'                                                            TUNIT5  = \'s       \'                                                            TTYPE6  = \'TSYS    \'           / system temperature in Kelvin                   TFORM6  = \'1D      \'                                                            TUNIT6  = \'K       \'                                                            COMMENT End of SDFITS CORE keywords/columns.                                    COMMENT Start of SDFITS DATA column and descriptive axes.                       TTYPE7  = \'DATA    \'           / actual data                                    TFORM7  = \'32768E  \'                                                            TUNIT7  = \'        \'                                                            TTYPE8  = \'TDIM7   \'           / data dimensions of the array                   TFORM8  = \'16A     \'                                                            TUNIT8  = \'        \'                                                            TTYPE9  = \'TUNIT7  \'                                                            TFORM9  = \'6A      \'                                                            TUNIT9  = \'        \'                                                            TTYPE10 = \'CTYPE1  \'           / first data axis is frequency-like              TFORM10 = \'8A      \'                                                            TUNIT10 = \'Hz      \'                                                            TTYPE11 = \'CRVAL1  \'                                                            TFORM11 = \'1D      \'                                                            TUNIT11 = \'Hz      \'                                                            TTYPE12 = \'CRPIX1  \'                                                            TFORM12 = \'1D      \'                                                            TUNIT12 = \'        \'                                                            TTYPE13 = \'CDELT1  \'                                                            TFORM13 = \'1D      \'                                                            TUNIT13 = \'Hz      \'                                                            TTYPE14 = \'CTYPE2  \'           / second axis is longitude-like                  TFORM14 = \'4A      \'                                                            TUNIT14 = \'        \'                                                            TTYPE15 = \'CRVAL2  \'                                                            TFORM15 = \'1D      \'                                                            TUNIT15 = \'deg     \'                                                            TTYPE16 = \'CTYPE3  \'           / third axis is latitude-like                    TFORM16 = \'4A      \'                                                            TUNIT16 = \'        \'                                                            TTYPE17 = \'CRVAL3  \'                                                            TFORM17 = \'1D      \'                                                            TUNIT17 = \'deg     \'                                                            CTYPE4  = \'STOKES  \'           / fourth axis is Stokes                          TTYPE18 = \'CRVAL4  \'                                                            TFORM18 = \'1I      \'                                                            TUNIT18 = \'        \'                                                            COMMENT End of SDFITS DATA column and descriptive axes.                         COMMENT Start of SDFITS SHARED keywords/columns.                                TTYPE19 = \'OBSERVER\'           / name of observer(s)                            TFORM19 = \'32A     \'                                                            TUNIT19 = \'        \'                                                            TTYPE20 = \'OBSID   \'           / observation description                        TFORM20 = \'32A     \'                                                            TUNIT20 = \'        \'                                                            PROJID  = \'TGBT17A_506_11\'     / project identifier                             TTYPE21 = \'SCAN    \'           / scan number                                    TFORM21 = \'1J      \'                                                            TUNIT21 = \'        \'                                                            TTYPE22 = \'OBSMODE \'           / observing mode                                 TFORM22 = \'32A     \'                                                            TUNIT22 = \'        \'                                                            TTYPE23 = \'FRONTEND\'           / frontend device                                TFORM23 = \'16A     \'                                                            TUNIT23 = \'        \'                                                            BACKEND = \'VEGAS   \'           / backend device                                 TTYPE24 = \'TCAL    \'           / calibration temperature                        TFORM24 = \'1E      \'                                                            TUNIT24 = \'K       \'                                                            TTYPE25 = \'VELDEF  \'           / velocity definition and frame                  TFORM25 = \'8A      \'                                                            TUNIT25 = \'        \'                                                            TTYPE26 = \'VFRAME  \'           / radial velocity of the reference frame         TFORM26 = \'1D      \'                                                            TUNIT26 = \'m/s     \'                                                            TTYPE27 = \'RVSYS   \'           / radial velocity, Vsource - Vtelescope          TFORM27 = \'1D      \'                                                            TUNIT27 = \'m/s     \'                                                            TTYPE28 = \'OBSFREQ \'           / observed center frequency                      TFORM28 = \'1D      \'                                                            TUNIT28 = \'Hz      \'                                                            TTYPE29 = \'LST     \'           / LST at midpoint of integration/scan            TFORM29 = \'1D      \'                                                            TUNIT29 = \'s       \'                                                            TTYPE30 = \'AZIMUTH \'           / azimuth                                        TFORM30 = \'1D      \'                                                            TUNIT30 = \'deg     \'                                                            TTYPE31 = \'ELEVATIO\'           / elevation                                      TFORM31 = \'1D      \'                                                            TUNIT31 = \'deg     \'                                                            TTYPE32 = \'TAMBIENT\'           / ambient temperature                            TFORM32 = \'1D      \'                                                            TUNIT32 = \'K       \'                                                            TTYPE33 = \'PRESSURE\'           / ambient pressure                               TFORM33 = \'1D      \'                                                            TUNIT33 = \'mmHg    \'                                                            TTYPE34 = \'HUMIDITY\'           / relative humidity                              TFORM34 = \'1D      \'                                                            TUNIT34 = \'        \'                                                            SITELONG=        -7.983983E+01 / E. longitude of intersection of the az/el axes SITELAT =         3.843312E+01 / N. latitude of intersection of the az/el axes  SITEELEV=         8.245950E+02 / height of the intersection of az/el axes       TTYPE35 = \'RESTFREQ\'           / rest frequency at band center                  TFORM35 = \'1D      \'                                                            TUNIT35 = \'Hz      \'                                                            TTYPE36 = \'FREQRES \'           / frequency resolution                           TFORM36 = \'1D      \'                                                            TUNIT36 = \'Hz      \'                                                            COMMENT End of SDFITS SHARED keywords/columns.                                  COMMENT Start of GBT-specific keywords/columns.                                 TTYPE37 = \'EQUINOX \'           / equinox of selected coordinate reference frame TFORM37 = \'1D      \'                                                            TUNIT37 = \'        \'                                                            TTYPE38 = \'RADESYS \'           / Equitorial coordinate system name              TFORM38 = \'8A      \'                                                            TTYPE39 = \'TRGTLONG\'           / target longitude in coord. ref. frame          TFORM39 = \'1D      \'                                                            TUNIT39 = \'deg     \'                                                            TTYPE40 = \'TRGTLAT \'           / target latitude in coord. ref. frame           TFORM40 = \'1D      \'                                                            TUNIT40 = \'deg     \'                                                            TTYPE41 = \'SAMPLER \'           / sampler description (e.g., "A1" or "A1xA3"     TFORM41 = \'12A     \'                                                            TUNIT41 = \'        \'                                                            TTYPE42 = \'FEED    \'           / (signal) feed number                           TFORM42 = \'1I      \'                                                            TUNIT42 = \'        \'                                                            TTYPE43 = \'SRFEED  \'           / reference feed number                          TFORM43 = \'1I      \'                                                            TUNIT43 = \'        \'                                                            COMMENT Feed offsets ARE included in the CRVAL2 and CRVAL3 columns              TTYPE44 = \'FEEDXOFF\'           / feed XEL offset                                TFORM44 = \'1D      \'                                                            TUNIT44 = \'deg     \'                                                            TTYPE45 = \'FEEDEOFF\'           / feed EL offset                                 TFORM45 = \'1D      \'                                                            TUNIT45 = \'deg     \'                                                            TTYPE46 = \'SUBREF_STATE\'       / subreflector state (1,0,-1) - 0=moving         TFORM46 = \'1I      \'                                                            TUNIT46 = \'        \'                                                            TTYPE47 = \'SIDEBAND\'           / resulting sideband (\'U\'pper or \'L\'ower)        TFORM47 = \'1A      \'                                                            TUNIT47 = \'        \'                                                            TTYPE48 = \'PROCSEQN\'           / scan sequence number                           TFORM48 = \'1I      \'                                                            TUNIT48 = \'        \'                                                            TTYPE49 = \'PROCSIZE\'           / number of scans in procedure                   TFORM49 = \'1I      \'                                                            TUNIT49 = \'        \'                                                            TTYPE50 = \'PROCSCAN\'           / scan\'s role in the procedure                   TFORM50 = \'16A     \'                                                            TUNIT50 = \'        \'                                                            TTYPE51 = \'PROCTYPE\'           / type of the procedure                          TFORM51 = \'16A     \'                                                            TUNIT51 = \'        \'                                                            TTYPE52 = \'LASTON  \'           / last \'on\' for position switching               TFORM52 = \'1J      \'                                                            TUNIT52 = \'        \'                                                            TTYPE53 = \'LASTOFF \'           / last \'off\' for position switching              TFORM53 = \'1J      \'                                                            TUNIT53 = \'        \'                                                            TTYPE54 = \'TIMESTAMP\'          / date and time of scan start                    TFORM54 = \'22A     \'                                                            TUNIT54 = \'UTC     \'                                                            TTYPE55 = \'QD_XEL  \'           / QuadrantDetector cross-elevation offset        TFORM55 = \'1D      \'                                                            TUNIT55 = \'deg     \'                                                            TTYPE56 = \'QD_EL   \'           / QuadrantDetector elevation offset              TFORM56 = \'1D      \'                                                            TUNIT56 = \'deg     \'                                                            TTYPE57 = \'QD_BAD  \'           / QuadrantDetector flag: 0=good,1=bad            TFORM57 = \'1I      \'                                                            TUNIT57 = \'        \'                                                            TTYPE58 = \'QD_METHOD\'          / Quad. Det. method A,B,C. Blank indicates none. TFORM58 = \'1A      \'                                                            TUNIT58 = \'        \'                                                            TTYPE59 = \'VELOCITY\'           / line velocity in rest frame                    TFORM59 = \'1D      \'                                                            TUNIT59 = \'m/s     \'                                                            TTYPE60 = \'ZEROCHAN\'           / zero channel                                   TFORM60 = \'1E      \'                                                            TUNIT60 = \'        \'                                                            TTYPE61 = \'DOPFREQ \'           / Doppler tracked frequency                      TFORM61 = \'1D      \'                                                            TUNIT61 = \'Hz      \'                                                            TTYPE62 = \'ADCSAMPF\'           / VEGAS ADC sampler frequency                    TFORM62 = \'1D      \'                                                            TUNIT62 = \'        \'           /                                                COMMENT Columns describing VEGAS spur locations                                 COMMENT   SPUR_CHANNEL = (J-VSPRVAL)*VSPDELT+VSPRPIX                            COMMENT   0 <= J <= 32                                                          COMMENT   spur channels are along the frequency axis                            TTYPE63 = \'VSPDELT \'           / channel increment between adjacent VEGAS spurs TFORM63 = \'1D      \'                                                            TUNIT63 = \'        \'                                                            TTYPE64 = \'VSPRVAL \'           / VEGAS spur number at VSPRPIX                   TFORM64 = \'1D      \'                                                            TUNIT64 = \'        \'                                                            TTYPE65 = \'VSPRPIX \'           / channel number of VEGAS spur VSPRVAL           TFORM65 = \'1D      \'                                                            TUNIT65 = \'        \'                                                            TTYPE66 = \'SIG     \'           / signal is true, reference is false             TFORM66 = \'1A      \'                                                            TUNIT66 = \'        \'                                                            TTYPE67 = \'CAL     \'           / cal ON is true, cal OFF is false               TFORM67 = \'1A      \'                                                            TUNIT67 = \'        \'                                                            TTYPE68 = \'CALTYPE \'           / LOW or HIGH, may eventually be other types     TFORM68 = \'8A      \'                                                            TTYPE69 = \'TWARM   \'           / 4mm RX ambient load temp (K)                   TFORM69 = \'1E      \'                                                            TUNIT69 = \'K       \'                                                            TTYPE70 = \'TCOLD   \'           / 4mm RX cold load temp (K)                      TFORM70 = \'1E      \'                                                            TUNIT70 = \'K       \'                                                            TTYPE71 = \'CALPOSITION\'        / 4mm RX table position                          TFORM71 = \'16A     \'                                                            TTYPE72 = \'IFNUM   \'           / Spectral window (IF) number                    TFORM72 = \'1I      \'                                                            TTYPE73 = \'PLNUM   \'           / Polarization number                            TFORM73 = \'1I      \'                                                            TTYPE74 = \'FDNUM   \'           / Feed number                                    TFORM74 = \'1I      \'                                                            COMMENT End of GBT-specific keywords/columns.                                   EXTNAME = \'SINGLE DISH\'        / name of this binary table extension            END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             '
    
    header = fits.Header.fromstring(table_header)
    
    return header


def read_sdfits(filename, ext='SINGLE DISH'):
    """
    """
    
    table, head, phead = _read_sdfits(filename, ext=ext)
    sdfits = SDFITS(table, head, phead=phead)
    
    return sdfits


def _read_sdfits(filename, ext='SINGLE DISH'):
    """
    """
    
    # Open the fits file.
    hdu = fits.open(filename)
    # Save the primary HDU header.
    phead = hdu[0].header
    
    # Find how many extensions match `ext`.
    sdhdus = []
    for i in range(len(hdu)):
        if hdu[i].name == ext:
            sdhdus.append(i)
    
    # Save the matching headers and tables.
    nhdu = len(sdhdus)
    tables = np.empty(nhdu, dtype=object)
    heads = np.empty(nhdu, dtype=object)
    for i,sdh in enumerate(sdhdus):
        tables[i] = hdu[sdh].data
        heads[i] = hdu[sdh].header           
    
    return tables, heads, phead
    

def write_sdfits(filename, sdfits, overwrite=False):
    """
    """
    
    table = sdfits.table
    header = sdfits.header
    
    phead = sdfits.phead
    empty_primary = fits.PrimaryHDU(header=phead)
    
    table_hdus = []
    for i,tab in enumerate(table):
        table_hdus.append(fits.BinTableHDU(data=tab, header=header[i], name='SINGLE DISH'))
    
    new_hdul = fits.HDUList([empty_primary,] + table_hdus)
    new_hdul.writeto(filename, overwrite=overwrite)
    #_write_sdfits(filename, table, header, overwrite=overwrite)
    
    
def _write_sdfits(filename, table, header, overwrite=False):
    """
    """
    
    new_hdu = fits.BinTableHDU.from_columns(table)
    
    for k,v in header.items(): 
        if k not in new_hdu.header.keys(): 
            new_hdu.header[k] = v
            
    new_hdu.writeto(filename, overwrite=overwrite)
    

def write_new_sdfits(filename, table, overwrite=False):
    """
    """
    
    header = make_sdfits_primary_header(date=None)
    empty_primary = fits.PrimaryHDU(header=header)
    
    header = make_sdfits_table_header()
    table_hdu = fits.BinTableHDU(data=table, header=header, name='SINGLE DISH')
    
    new_hdul = fits.HDUList([empty_primary, table_hdu])
    new_hdul.writeto(filename, overwrite=overwrite)
    
    
