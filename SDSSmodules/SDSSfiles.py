# This module contains functions, which read files, txt as well as FITS

import numpy as np

# astropy is used to work with fits files
from astropy.io import fits

# contains astropy coordinate objects used to work with coordinates
from astropy.coordinates import SkyCoord
from astropy import units as u

# fitsio based on c
import fitsio
from fitsio import FITS, FITSHDR

# re contains the function search, which is used to search for 
# spec and spSpec in the filenames
import re

from SDSSmodules.SDSScorrections import *

# import rebin_1d from congrid
from congrid import rebin_1d

#####################################################################################
################################### read files ######################################
#####################################################################################



def read_arrays_and_interpolate(files, wave_low, wave_high, npix = 4000, corr_flag = False):
    # files: list of FITS files
    # wave_low:  beginning wavelength from where to interpolate
    # wave_high: end wavelenght where to interpolate to
    # npix:      number of pixels to interpolate to (default 4000)
    # corr_flag: flag which if true, we also read the flux_correction
    #            function 
    
    # this function works in the folowing way:
    # 1) it receives a list of FITS files
    # 2) it checks the filetype of the FITS files
    # 3) it reads the FITS files in a loop using the corresponding read_*_fitsio
    #    function
    # 4) we obtain the indices from the wavelength array which we need to consider
    #    in the flux arrays to interpolate them
    #    We also check for isfinite in flux and flux_corr
    # 5) using rebin_1d we interpolate all arrays to exactly the same length
    # 6) after 10000 files, we dump the data to a cPickle file
    ########## WARNING: We immediately drop spectra, which have wavelength arrays with
    ################### np.min(wave) > wave_low 
    ################### np.min(wave) < wave_high

    # variable names:
    # flux, flux_corr: LISTS (already interpolated!!!)
    # fl, f_corr: arrays of ONE single spectrum
    

    filetype = check_filetype(files[0])
    
    for i, file in enumerate(files):
        # if we are not on the first iteration and i mod 5000 == 0
        # then we dump the lists to file
        if i is not 0 and i % 5000 == 0:
            # open file to dump data into
            outfile_name = '/mnt/Daten/Uni/SDSS/data/pickled_data/data_' + str(filetype) + '_' + str(i) + '.dat' 
            outfile = open(outfile_name, 'wb')
            cPickle.dump(flux,       outfile, -1)
            cPickle.dump(flux_error, outfile, -1)
            if corr_flag == True:
                cPickle.dump(flux_corr, outfile, -1)
            outfile.close()
            
            # now we reset the lists
            flux          = []
            flux_error    = []
            wave_min_list = []
            if corr_flag == True:
                flux_corr  = []
        elif i is 0:
            flux          = []
            flux_error    = []
            wave_min_list = []
            if corr_flag == True:
                flux_corr  = []

        if i % 500 == 0:
            print i, 'spectra read'

        spec = spectrum()        
        if filetype == 1:
            read_spSpec_fitsio(files[i], spec, None)
        if filetype == 2:
            read_spec_fitsio(files[i], spec, None)
        if filetype == 3:
            read_speclya_fitsio(files[i], spec, None)
        if filetype == 4:
            read_mockspec_fitsio(files[i], spec)
    
        wave   = spec.wave
        if corr_flag == True:
            # this is where it gets ugly? 
            # if corr_flag is true, we need to 
            # if perform_corr is set to 0, we simply read the flux correction from the corrected FITS files
            # we don't have a list of files in /home/z5010843/data/flux_calib_duplicates/
            # so we just try to open the current files
            # read correction function
            file_path = '/mnt/Daten/Uni/SDSS/data/DR12_flux_calib/' + files[i]
            try:
                hdu            = fitsio.FITS(file_path)
                spec.flux_corr = hdu[1]['corr'][:]
                hdu.close()
            except IOError:
                # set flag to 1, so that we skip this duplicate
                # if an IOError is called, we probably miss the file containing the flux correction
                # we could in principle then call perform_flux_correction_adaptive and get it
                # Do we want that?
                continue
            except ValueError:
                # set flag to 1, so that we skip this duplicate
                # if this happens, the file doesn't contain the corr column it seems? broken file?
                continue
            f_corr = spec.flux_corr / spec.flux
        if np.min(wave) > wave_low or np.max(wave) < wave_high:
            print('The wavelength array of %s is inside the wave_low and wave_high boundaries' % files[i])
            print 'wave_low  = ', wave_low
            print 'wave_high = ', wave_high
            print 'wave min and max:', np.min(wave), np.max(wave)
            print 'file: ', files[i]
            print 'Skip this file and continue'
            continue
       
        index = np.where( (wave > wave_low ) &
                          (wave < wave_high) &
                          (np.isfinite(spec.flux) == True) &
                          (np.isfinite(f_corr) == True))[0]
        # the last line should hopefully be redundant!
        # now we have the indices of all elements in the wavelength (and thus flux, flux_corr)
        # arrays, which we will now interpolate
        fl   = spec.flux[index]
        fl   = rebin_1d(fl,   npix, method='spline')
        wave = rebin_1d(wave, npix, method='spline')
        flux.append(fl)
        wave_min_list.append(wave)
        if corr_flag is True:
            f_corr = f_corr[index]
            f_corr = rebin_1d(f_corr, npix, method='spline')
            flux_corr.append(f_corr)
        
        # at this point we should be done for this spectrum
        ############### ITERATION END




def read_resid_corr(file):
    # This function reads the file containing the residual correction
    # for the Balmer problem in SDSS
    # returns an array containing the correction values for each
    # wavelength
    print 'start creating resid_corr array'

    resid_file = open(file, 'r').readlines()
    resid_file = resid_file[1:]
    resid_corr = []
    for line in resid_file:
        wave, resid = line.split()
        resid_corr.append(float(resid))

    resid_corr = np.asarray(resid_corr)
    print 'finished resid_corr array'
    return resid_corr


def check_filetype(filename):
# Function which checks the filetype of a given file (string)
# simply searches for spSpec and spec in the filename 
    if re.search("spSpec", filename):
        return 1
    elif re.search('mockspec', filename):
        return 4
    elif re.search("speclya", filename):
        return 3
    elif re.search("spec", filename):
        return 2
    else:
        return 0



################################### FITS files ######################################

def read_mockspec_fitsio(qso, spec):

    qso = qso.rstrip()
    hdu = fitsio.FITS(qso)

    # Coordinates:
    TableHDU1 = hdu[1]
    ra = TableHDU1['RA'][0]
    dec = TableHDU1['DEC'][0]
    spec.coordinates = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame='fk5')
    spec.MJD = TableHDU1['MJD'][0]
    spec.plateid = TableHDU1['plate'][0]
    spec.fiberid = TableHDU1['fiber'][0]
    # Table data from HDU2; contains redshift, MJD and psf magnitudes
    spec.z = TableHDU1['Z_VI'][0]
    zwarn  = TableHDU1['Z_WARNING'][0]


    # Can't find cpix in header, assume it's 1.
    spec.cpix = 1
    # defining hdu1 variable is faster than accessing hdu[1] each time
    hdu2header = hdu[2].read_header()
    spec.npix = hdu2header['NAXIS2']

    TableHDU2 = hdu[2]
    spec.flux = TableHDU2['flux'][:]
    spec.flux_error = TableHDU2['noise'][:]
    
    spec.model_sdss = TableHDU2['no_abs'][:]
    spec.wave    = 10.0**(TableHDU2['loglam'][:])

#    print spec.model_sdss



            
    del(TableHDU1)
    del(TableHDU2)
    hdu.close()

    # if zerr < 0 or zwarn != 0:
    #     if settings is not None and settings.debug == 1:
    #         print 'bad zerr'
    #         print zerr, zwarn, spec.z
    #     spec.z = 999
    #     return 2

    return 0
    


def read_speclya_fitsio(qso, spec, settings=None):

# TODO: we could check if Z_ERR is < 0, if so, it means that the 
# estimate for the redshift is probably incorrect

# The same function as read_spec(), but based on fitsio, which is much faster, since
# it's based on cfitsio
# Function whic is called to read spec FITS file

    #Open FITS file
    qso = qso.rstrip()
    hdu = fitsio.FITS(qso)
    spec.filename = qso

    # Read in all necessary information from Header of
    # defining hdu0header is faster than accessing hdu[0].header each time
    hdu0header = hdu[0].read_header()

    spec.plateid = hdu0header['PLATEID']
    spec.beginwl = hdu0header['COEFF0']
    spec.deltawl = hdu0header['COEFF1']

    # Coordinates:
    ra = hdu0header['PLUG_RA']
    dec = hdu0header['PLUG_DEC']
    spec.coordinates = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame='fk5')
    spec.fiberid = hdu0header['FIBERID']
    # Can't find cpix in header, assume it's 1.
    spec.cpix = 1
    # defining hdu1 variable is faster than accessing hdu[1] each time
    hdu1header = hdu[1].read_header()
    spec.npix = hdu1header['NAXIS2']
    spec.PCA_qual = hdu1header['CONTFLG']
#    print 'PCA quality:'
#    print spec.PCA_qual

    # Get the table data from HDU 1; it contains the flux
    # Change tabledata = hdu.data to tabledata = table.view(np.recarray)
    TableHDU1 = hdu[1]
    # TODO: maybe make copy?
    # correct for DLA
    spec.flux = TableHDU1['flux'][:]*TableHDU1['DLA_CORR'][:]#.copy()

    # Error is given in inverse variance. To get STD, we have to take sqrt(1/ivar)
    # TODO: Check if need to copy data?
    # correct ivar values
    spec.flux_error = np.sqrt(1/TableHDU1['ivar'][:])/TableHDU1['NOISE_CORR'][:]
    
    spec.model_sdss = TableHDU1['CONT'][:]
    spec.mask_comb = TableHDU1['MASK_COMB'][:]

#    print spec.model_sdss

    # Table data from HDU2; contains redshift, MJD and psf magnitudes
    spec.z = hdu[2]['Z'][0]
    zerr   = hdu[2]['Z_ERR'][0]
    zwarn  = hdu[2]['ZWARNING'][0]
    spec.MJD = hdu[2]['MJD'][0]

    #print "%.15f" % spec.z

    # spec.mag = TableHDU2['PSFMAG'].copy(

    #TODO: Read in primary target information?

    # run over all wavelengths (pixels) and set the wavelength array
    # as well as the status array and signal to noise ratio
    # wave, status and snr are all numpy arrays
    spec.wave    = 10.0**(spec.beginwl + (np.arange(spec.npix) + 1 - spec.cpix)*spec.deltawl)
    # use Numpy's logical operations, to determine the status array
    spec.status  = spec.flux_error > 0
    # determine snr by multiplying the arrays
    # np.seterr(all='raise')
    # try:
    #     spec.snr     = spec.flux / spec.flux_error
    # except FloatingPointError:
    #     for i in xrange(spec.npix):
    #         if spec.flux_error[i] != 0:
    #             spec.snr.append(spec.flux[i] / spec.flux_error[i])
    #         else:
    #             spec.snr.append(0)
            
    del(TableHDU1)
    hdu.close()

    if zerr < 0 or zwarn != 0:
        if settings is not None and settings.debug == 1:
            print 'bad zerr'
            print zerr, zwarn, spec.z
        spec.z = 999
        return 2

    return 0



def read_spec_fitsio(qso, spec, settings=None, resid_corr = None):

# The same function as read_spec(), but based on fitsio, which is much faster, since
# it's based on cfitsio
# Function whic is called to read spec FITS file

    #Open FITS file
    qso = qso.rstrip()
    hdu = fitsio.FITS(qso)
    spec.filename = qso

    # Read in all necessary information from Header of
    # defining hdu0header is faster than accessing hdu[0].header each time
    hdu0header = hdu[0].read_header()

    spec.MJD = hdu0header['MJD']
    spec.plateid = hdu0header['PLATEID']
    spec.beginwl = hdu0header['COEFF0']
    spec.deltawl = hdu0header['COEFF1']

    # Coordinates:
    ra = hdu0header['PLUG_RA']
    dec = hdu0header['PLUG_DEC']
    spec.coordinates = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame='fk5')
    spec.fiberid = hdu0header['FIBERID']
    # Can't find cpix in header, assume it's 1.
    spec.cpix = 1
    # defining hdu1 variable is faster than accessing hdu[1] each time
    hdu1header = hdu[1].read_header()
    spec.npix = hdu1header['NAXIS2']

    # Get the table data from HDU 1; it contains the flux
    # Change tabledata = hdu.data to tabledata = table.view(np.recarray)
    TableHDU1 = hdu[1]
    spec.flux = TableHDU1['flux'][:]
    # Error is given in inverse variance. To get STD, we have to take sqrt(1/ivar)
    spec.flux_error = np.sqrt(1/TableHDU1['ivar'][:])
    if resid_corr is not None:
        apply_resid_corr(spec, resid_corr)

    #spec.sky = TableHDU1['sky'][:]
    #    spec.model_sdss = TableHDU1['model'][:]

    # need to create better mask
    spec.mask_comb = TableHDU1['and_mask'][:]


    # Table data from HDU2; contains redshift, MJD and psf magnitudes
    spec.z = hdu[2]['Z'][0]
    #print 'ok,', spec.z
    zerr   = hdu[2]['Z_ERR'][0]
    zwarn  = hdu[2]['ZWARNING'][0]

    # spec.mag = TableHDU2['PSFMAG']

    # run over all wavelengths (pixels) and set the wavelength array
    # as well as the status array and signal to noise ratio
    # wave, status and snr are all numpy arrays
    spec.wave    = 10.0**(spec.beginwl + (np.arange(spec.npix) + 1 - spec.cpix)*spec.deltawl)
    # use Numpy's logical operations, to determine the status array
    spec.status  = spec.flux_error > 0
    # determine snr by multiplying the arrays
    # np.seterr(all='raise')
    # try:
    #     spec.snr     = spec.flux / spec.flux_error
    # except FloatingPointError:
    #     for i in xrange(spec.npix):
    #         if spec.flux_error[i] != 0:
    #             spec.snr.append(spec.flux[i] / spec.flux_error[i])
    #         else:
    #             spec.snr.append(0)
            
    del(TableHDU1)
    hdu.close()

    if zerr < 0 or zwarn != 0:
        if settings is not None and settings.debug == 1:
            print 'bad zerr'
            print zerr, zwarn, spec.z
        spec.z = 999
        return 2

    return 0



def read_spSpec_fitsio(qso, spec, settings=None):
# same function as read_spSpec implemented with fitsio instead of astropy
# Function which is called to read spSpec FITS file

    #Open FITS file
    qso = qso.rstrip()
    hdu = fitsio.FITS(qso)

    #Read in all necessary information from Header of HDU 0:
    hdu0header = hdu[0].read_header()
    spec.plateid = hdu0header['PLATEID']
    spec.fiberid = hdu0header['FIBERID']

    # Coordinates:
    ra = hdu0header['RAOBJ']
    dec = hdu0header['DECOBJ']
    spec.coordinates = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame='fk5')
    spec.MJD = hdu0header['MJD']
    spec.z = hdu0header['Z']
    zerr   = hdu0header['Z_ERR']
    spec.mag = hdu0header['MAG']
    spec.npix = hdu0header['NAXIS1']
    spec.beginwl = hdu0header['CRVAL1']
    spec.deltawl = hdu0header['CD1_1']
    spec.cpix = hdu0header['CRPIX1']
    #Open ImageData from HDU 0 and thus retrieve the flux and the error
    ImageData = hdu[0].read()
    spec.flux = ImageData[0, 0:spec.npix].copy()
    spec.flux_error = ImageData[2, 0:spec.npix].copy()
    
    # run over all wavelengths (pixels) and set the wavelength array
    # as well as the status array and signal to noise ratio
    # wave, status and snr are all numpy arrays
    spec.wave    = 10.0**(spec.beginwl + (np.arange(spec.npix) + 1 - spec.cpix)*spec.deltawl)
    # use Numpy's logical operations, to determine the status array
    spec.status  = spec.flux_error > 0
    # determine snr by multiplying the arrays
    # np.seterr(all='raise')
    # try:
    #     spec.snr     = spec.flux / spec.flux_error
    # except FloatingPointError:
    #     for i in xrange(spec.npix):
    #         if spec.flux_error[i] != 0:
    #             spec.snr.append(spec.flux[i] / spec.flux_error[i])
    #         else:
    #             spec.snr.append(0)
    del(ImageData)
    hdu.close()

    if zerr < 0:
        print zerr, spec.z
        spec.z = 999
        return 2

    return 0
# What is primary target flag needed for?

def read_spec(qso, spec, settings = None):
# Function whic is called to read spec FITS file

    #Open FITS file
    hdu = fits.open(qso)
    # Read in all necessary information from Header of
    # defining hdu0header is faster than accessing hdu[0].header each time
    hdu0header = hdu[0].header

    spec.plateid = hdu0header['PLATEID']
    spec.fiberid = hdu0header['FIBERID']
    spec.MJD = hdu0header['MJD']

    # Coordinates:
    ra = hdu0header['PLUG_RA']
    dec = hdu0header['PLUG_DEC']
    spec.coordinates = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame='fk5')
    # Check if we make cuts on coordinates. If so, check if object lies within coordinate window
    # if settings.coords:
    #     if (spec.coordinates.galactic.l.deg >= settings.l_min and
    #         spec.coordinates.galactic.l.deg <= settings.l_max and
    #         spec.coordinates.galactic.b.deg >= settings.b_min and
    #         spec.coordinates.galactic.b.deg <= settings.b_max) == 0:
    #         # return 1 if we neglect this quasar
    #         return 1
    
    spec.beginwl = hdu0header['COEFF0']
    spec.deltawl = hdu0header['COEFF1']
    # Can't find cpix in header, assume it's 1.
    spec.cpix = 1
    # defining hdu1 variable is faster than accessing hdu[1] each time
    hdu1 = hdu[1]
    spec.npix = hdu1.header['NAXIS2']
    # Get the table data from HDU 1; it contains the flux
    # Change tabledata = hdu.data to tabledata = table.view(np.recarray)
    TableHDU1 = hdu1.data
    spec.flux = TableHDU1['flux'].copy()

    # Error is given in inverse variance. To get STD, we have to take sqrt(1/ivar)
    # TODO: Check if need to copy data?
    spec.flux_error = np.sqrt(1/TableHDU1['ivar'])

    # Table data from HDU2; contains redshift, MJD and psf magnitudes
    spec.z = hdu[2].data['Z'][0]

    # spec.mag = TableHDU2['PSFMAG'].copy(

    #TODO: Read in primary target information?

    # run over all wavelengths (pixels) and set the wavelength array
    # as well as the status array and signal to noise ratio
    # wave, status and snr are all numpy arrays
    spec.wave    = 10.0**(spec.beginwl + (np.arange(spec.npix) + 1 - spec.cpix)*spec.deltawl)
    # use Numpy's logical operations, to determine the status array
    spec.status  = spec.flux_error > 0
    # determine snr by multiplying the arrays
    # np.seterr(all='raise')
    # try:
    #     spec.snr     = spec.flux / spec.flux_error
    # except FloatingPointError:
    #     for i in xrange(spec.npix):
    #         if spec.flux_error[i] != 0:
    #             spec.snr.append(spec.flux[i] / spec.flux_error[i])
    #         else:
    #             spec.snr.append(0)

    del(TableHDU1)
    hdu.close()
    return 0


def read_spSpec(qso, spec, settings=None):
# Function which is called to read spSpec FITS file

    #Open FITS file
    hdu = fits.open(qso)
    #Read in all necessary information from Header of HDU 0:
    hdu0header = hdu[0].header
    spec.plateid = hdu0header['PLATEID']
    spec.fiberid = hdu0header['FIBERID']

    # Coordinates:
    ra = hdu0header['RAOBJ']
    dec = hdu0header['DECOBJ']
    spec.coordinates = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame='fk5')
    # Check if we make cuts on coordinates. If so, check if object lies within coordinate window
    # if settings.coords:
    #     if (spec.coordinates.galactic.l.deg >= settings.l_min and
    #         spec.coordinates.galactic.l.deg <= settings.l_max and
    #         spec.coordinates.galactic.b.deg >= settings.b_min and
    #         spec.coordinates.galactic.b.deg <= settings.b_max) == 0:
    #         # return 1 if we neglect this quasar
    #         print spec.coordinates.galactic.l.deg, spec.coordinates.galactic.b.deg
    #         print settings.l_min, settings.b_min
    #         return 1


    spec.MJD = hdu0header['MJD']
    spec.z = hdu0header['Z']
    spec.mag = hdu0header['MAG']
    spec.npix = hdu0header['NAXIS1']
    spec.beginwl = hdu0header['CRVAL1']
    spec.deltawl = hdu0header['CD1_1']
    spec.cpix = hdu0header['CRPIX1']
    #Open ImageData from HDU 0 and thus retrieve the flux and the error
    ImageData = hdu[0].data.copy()
    spec.flux = ImageData[0, 0:spec.npix].copy()
    # 2 instead of ImageData[1,..
    spec.flux_error = ImageData[2, 0:spec.npix].copy()
    
    # run over all wavelengths (pixels) and set the wavelength array
    # as well as the status array and signal to noise ratio
    # wave, status and snr are all numpy arrays
    spec.wave    = 10.0**(spec.beginwl + (np.arange(spec.npix) + 1 - spec.cpix)*spec.deltawl)
    # use Numpy's logical operations, to determine the status array
    spec.status  = spec.flux_error > 0
    # determine snr by multiplying the arrays
    # np.seterr(all='raise')
    # try:
    #     spec.snr     = spec.flux / spec.flux_error
    # except FloatingPointError:
    #     for i in xrange(spec.npix):
    #         if spec.flux_error[i] != 0:
    #             spec.snr.append(spec.flux[i] / spec.flux_error[i])
    #         else:
    #             spec.snr.append(0)

    del(ImageData)
    hdu.close()
    return 0
# What is primary target flag needed for?


################################### additional information from FITS files ##########

def get_array_from_ind_exposures(spec, settings, return_array = 0, check_weather = 0):

    hdu = fitsio.FITS(spec.filename)
    up  = hdu[-1].get_extnum() + 1
    # pressure, temperature, humidity
    header = hdu[0].read_header()
    #P = float(header['pressure'])
    #T = float(header['AIRTEMP'])
    #h = float(header['HUMIDITY'])

    altitude = []
    azimuth  = []
    seeing50 = []
    seeing80 = []
    pressure = []
    airtemp  = []
    humidity = []
    snr      = []
    flux     = []
    error    = []
    wave     = []
    spec.temp_flag = 0
    
    for i in xrange(4, up):
        header = hdu[i].read_header()
        seeing50.append(float(header['SEEING50']))
        seeing80.append(float(header['SEEING80']))
        if check_weather == 1:
            try:
                pressure.append(float(header['pressure']))
            except ValueError:
                raise ValueError('No weather data found in file!')
        else:
            pressure.append(float(header['pressure']))
        azimuth.append(float(header['AZ']))
        altitude.append(float(header['ALT']))
        humidity.append(float(header['humidity']))
        try:
            airtemp.append(float(header['airtemp']))
        except ValueError:
            try:
                airtemp.append(float(header['temp']))
                if np.isnan(airtemp[i-4]) == True:
                    # in this case we set flag and assume 17.5 Celsius
                    #print 'temp', spec.file
                    spec.temp_flag = 1
                    spec.T = 17.5
            except ValueError:
                print 'not called temp either. not available?'
                import sys
                sys.exit()
        wave.append(np.around(hdu[i]['loglam'][:], decimals=4))
        flux.append(hdu[i]['flux'][:])
        error.append(1/hdu[i]['ivar'][:])
        snr.append(np.mean(flux[-1]/error[-1]))

    seeing50_array = np.asarray(seeing50)
    seeing80_array = np.asarray(seeing80)
    altitude_array = np.asarray(altitude)
    azimuth_array  = np.asarray(altitude)
    pressure_array = np.asarray(pressure)
    airtemp_array  = np.asarray(airtemp)
    humidity_array = np.asarray(humidity)
    snr_array      = np.asarray(snr)
    
    loglam = hdu[1]['loglam'][:]

    seeing50 = np.mean(seeing50_array)
    seeing80 = np.mean(seeing80_array)
    altitude = np.mean(altitude_array)
    azimuth = np.mean(azimuth_array)
    P = np.mean(pressure_array)
    if spec.temp_flag == 0:
        T = np.mean(airtemp_array)
    else:
        T = spec.T

    h = np.mean(humidity_array)
    # Convert pressure from Inches of Hg to HPa
    P = P * 3386.389 / 100
    try:
        ppW = partial_water_pressure(P, T, h)
    except ValueError:
        print 'Partial water pressure not successfully calculated'
        raise
    if np.isnan(ppW):
        print 'ppW', spec.filename

    hdu.close()

    if return_array == 0:
        spec.alt = altitude
        spec.az  = azimuth
        spec.P   = P
        if spec.temp_flag == 0:
            spec.T   = T
        spec.h   = h
        spec.ppW = ppW
        return seeing50, seeing80, altitude, P, T, ppW
    elif return_array == 1:
        return seeing50_array, seeing80_array, altitude_array, pressure_array, airtemp_array, humidity_array, snr_array, flux, error, wave, loglam

