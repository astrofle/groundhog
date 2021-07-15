===============
Groundhog tests
===============

GBTIDL tests
============

Many of the `Groundhog` functions are tested against the output of `GBTIDL`.
Here we list the commands used to generate the test files.

AGBT19B_332_04_scan5_ifnum4_intnum1_plnum0.ascii
------------------------------------------------

In GBTIDL:
   * filein,'AGBT19B_334_04_3C353.raw.vegas.A.fits'
   * getps,5,ifnum=4,intnum=1,plnum=0
   * Then, using the plotter window, save the spectra as ascii.
   

TGBT20A_506_01_scan6_ifnum0_plnum0.ascii
----------------------------------------

In GBTIDL:
   * offline,'TGBT20A_506_01'
   * getps,6,ifnum=0,plnum=0
   * Then, using the plotter window, save the spectra as ascii.
   

NGC6946_scan_14_26_window0_feed0_pol0.fits
------------------------------------------

On a GBO machine:
   * gbtpipeline -i /home/scratch/dfrayer/DATAdemo/TGBT17A_506_11.raw.vegas -m "14:26" --refscan 27
