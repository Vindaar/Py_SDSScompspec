# This module contains the functions, which write output to disk

import numpy as np

# contains function to call date and time
from datetime import datetime
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

################################################################################
############################## Output ##########################################
################################################################################

def build_fits_file(cspec, spec, outfile, settings):
# Build the actual FITS output file
    # set a few variables
    i = 0
    coeff1 = 0.0001
    nspec = len(spec)
    # TODO: FOLLOWING line is wrong? Should be nhist[i]. Changed!
    # while cspec.nhist == 0:
    while cspec.nhist[i] == 0:
        i += 1
    sidx = i
    npix = 0
    date_time = datetime.now()

    hdu0_row1 = []              # 1 - average transmittance metal calibrated
    hdu0_row2 = []              # error of average transmittance
    hdu0_row3 = []              # number of contributing objects to wavelength
    hdu0_row4 = []              # 1 - average transmittance not calibrated
    # changed to 4000, since sidx pretty big in DR10
    for i in xrange(sidx, 4000):
        if cspec.nhist[i] > 0:
            npix = i+1


    # Add a few header keys
    prihdr = fits.Header()
    prihdr['PROGRAM'] = settings.program_name
    prihdr['AUTHOR']  = ('Sebastian Schmidt', 'derived from SDSScompspec C code by Jon Oullet')
    prihdr['DATE']    = (date_time.strftime("%A, %d. %B %Y %I:%M%p"), 'Date created')
    prihdr['ARRAY0']  = ('DA', '1 - Flux / Flux_count')
    prihdr['ARRAY1']  = ('ERR', '1 Sigma')
    prihdr['ARRAY2']  = ('NPOINTS', 'Number of contributing pixels')
    prihdr['NSPECTRA']= (cspec.spectra_count, 'Total number of used / contributing QSOs')
    prihdr['NSPEC']   = (nspec, 'Total number of QSOs considered')
    prihdr['VACUUM']  = (1, 'Wavelengths are in vacuum')
    prihdr['DC_FLAG'] = (1, 'Log-Linear Flag')
    prihdr['CRPIX1']  = (1, 'Starting pixel (1-indexed)')
    prihdr['COEFF0']  = np.log10(cspec.wave[sidx])
    prihdr['COEFF1']  = coeff1
    prihdr['CRVAL1']  = np.log10(cspec.wave[sidx])
    prihdr['CD1_1']   = coeff1
    prihdr['DR']      = (10, 'SDSS Data release used')

    # run over all pixels and append elements to new arrays
    # note: row1 is technically not the 'optical depth tau'
    # instead it is 1 - average transmittance
    # row1 will be calibrated for metal afterwards
    # row2 simply the error of the average transmittance
    # row3 the number of contributing objects to each wavelength
    # row4 is the average transmittance not calibrated for metals
    print 'creation of arrays for HDU0'
    for i in xrange(npix):
        hdu0_row1.append(1 - cspec.flux[sidx+i])
        hdu0_row2.append(cspec.flux_error[sidx+i])
        hdu0_row3.append(cspec.nhist[sidx+i])
    # create row4 of the hdu. At the moment the same as row1, however row1 
    # will be changed after this
    hdu0_row4 = list(hdu0_row1)

    # Next up: metal corrections, based on:
    # \tau_corr(z) = 0.72(1+z)^(0.17) \tau_uncorr(z)
    # Schaye et al. (2003)
    # we need to assign a redshift to each pixel, thus we create the
    # zem array. Then, we first calculate the real optical depth from
    # hdu0_row1 (which technically is 1 - average transmittance) by
    # tau = - log10(average transmittance)
    # apply the correction and revert the optical depth back to 
    # 1 - average transmittance
    np.seterr(all='raise')
    print 'applying metal correction'
    zem = 10**(np.log10(cspec.wave[sidx]) + coeff1*np.arange(npix)) / 1215.6701 - 1
    # TODO: change len(hdu0_row1) to npix? 
    for i in xrange(len(hdu0_row1)):
        try:
            hdu0_row1[i] = -np.log10(1 - hdu0_row1[i])
            zem_temp     = zem[i]
            hdu0_row1[i] = 0.72*(1+zem_temp)**(0.17) * hdu0_row1[i]
            hdu0_row1[i] = 1 - 10**(-hdu0_row1[i])
        except FloatingPointError:
            print i, 'set to zero'
            # TODO: probably shouldn't be zero? 1 instead, because that's when the error occurs
            hdu0_row1[i] = 1

    # zip arrays; will create one array of 3 tuples from the three
    # individual arrays
    print 'zipping arrays for HDU0'
    zipped_hdu0 = zip(hdu0_row1, hdu0_row2, hdu0_row3, hdu0_row4)
    # currently wrong format, so we transpose the array
    zipped_hdu0 = np.transpose(zipped_hdu0)

    # create the Primary ouput HDU and print it to file
    hdu = fits.PrimaryHDU(zipped_hdu0, header=prihdr)

    # Create 2nd HDU containing information on all QSOs
    # Calculate arrays, which will fill the table HDU
    # TODO: optimize creation of arrays! Lots of potential probably

    print 'creation of arrays for HDU1'
    mjd_array         = np.zeros(nspec)
    plate_array       = np.zeros(nspec)
    fiber_array       = np.zeros(nspec)
    alpha_array       = np.zeros(nspec)
    alpha_error_array = np.zeros(nspec)
    chisq_array       = np.zeros(nspec)
    ra_array          = np.zeros(nspec)
    dec_array         = np.zeros(nspec)
    l_array           = np.zeros(nspec)
    b_array           = np.zeros(nspec)
    zem_array         = np.zeros(nspec)
    flag_array        = np.zeros(nspec)
    smag_array        = []

    for i in xrange(nspec):
        mjd_array[i]         = spec[i].MJD
        plate_array[i]       = spec[i].plateid
        fiber_array[i]       = spec[i].fiberid
        alpha_array[i]       = spec[i].alpha
        alpha_error_array[i] = spec[i].alpha_error
        chisq_array[i]       = spec[i].chisq
        ra_array[i]          = spec[i].coordinates.ra.degree
        dec_array[i]         = spec[i].coordinates.dec.degree
#        l_array[i]           = spec[i].coordinates.galactic.l.degree
#        b_array[i]           = spec[i].coordinates.galactic.b.degree
        zem_array[i]         = spec[i].z
        flag_array[i]        = int(spec[i].flag)
        smag_array.append(spec[i].smag)
        if i % 500 == 0:
            print i, "spectra added to arrays"
    
    print "create galactic coordinate arrays"
    coord_objs = SkyCoord(ra = ra_array*u.degree, dec = dec_array*u.degree, frame='fk5')
    l_array    = coord_objs.galactic.l.deg
    b_array    = coord_objs.galactic.b.deg

    # mjd_array         = map(lambda spec: int(spec.MJD), spec)
    # plate_array       = map(lambda spec: int(spec.plateid), spec)
    # fiber_array       = map(lambda spec: int(spec.fiberid), spec) 
    # alpha_array       = map(lambda spec: spec.alpha, spec)
    # alpha_error_array = map(lambda spec: spec.alpha_error, spec)
    # ra_array          = map(lambda spec: spec.coordinates.ra.degree, spec)
    # dec_array         = map(lambda spec: spec.coordinates.dec.degree, spec)
    # l_array           = map(lambda spec: spec.coordinates.galactic.l.degree, spec)
    # b_array           = map(lambda spec: spec.coordinates.galactic.b.degree, spec)
    # zem_array         = map(lambda spec: spec.z, spec)
    # smag_array        = map(lambda spec: spec.smag, spec)
    # flag_array        = map(lambda spec: int(spec.flag), spec)

    # write arrays to Table HDU:
    print 'write arrays to HDU1'
    TableHDU = fits.BinTableHDU.from_columns(
        fits.ColDefs([fits.Column(name='MJD',    format='J', array=mjd_array),
                      fits.Column(name='PLATE',  format='J', array=plate_array),
                      fits.Column(name='FIBER',  format='J', array=fiber_array),
                      fits.Column(name='RA',     format='D', array=ra_array),
                      fits.Column(name='DEC',    format='D', array=dec_array),
                      fits.Column(name='ALPHA',  format='D', array=alpha_array),
                      fits.Column(name='EALPHA', format='D', array=alpha_error_array),
                      fits.Column(name='ChiSq',  format='D', array=chisq_array),
                      fits.Column(name='l',      format='D', array=l_array),
                      fits.Column(name='b',      format='D', array=b_array),
                      fits.Column(name='ZEM',    format='D', array=zem_array),
                      fits.Column(name='MAGSPEC',format='PD()', array=smag_array),
                      fits.Column(name='FLAG',   format='I', array=flag_array)]))
    # TODO: check if I can change smag column a little bit


    # Now write new keywords to Table header
    table_header = TableHDU.header
    table_header['NSPEC']  = (nspec, 'Total number of Quasars')
    table_header['MEAN_A'] = (cspec.mean_a, 'Mean spectral index')
    table_header['SIGMA_A']= (cspec.sigma_a, 'Standard deviation on alpha')
    table_header['MED_A']  = (cspec.median_a, 'Median spectral index')
    table_header['SIQR_A'] = (cspec.siqr_a, '68% semi-interquartile range on alpha')

    

    # hdu.writeto('test.fits')
    print 'write FITS file'
    hdulist = fits.HDUList([hdu, TableHDU])
    # write to outfile and 'clobber=True' -> overwrite if existing
    hdulist.writeto(outfile, clobber=True)



def write_compspectra_to_file(compspectra, settings, outfile):
    # TODO: need extra stuff, but for now should work

# Build the actual FITS output file
    # set a few variables
    i = 0
    coeff1 = 0.0001
    cspec     = compspectra[0]
    ncspec    = np.size(compspectra)
    npix      = cspec.npix
    date_time = datetime.now()

    # Add a few header keys
    prihdr = fits.Header()
    prihdr['PROGRAM'] = settings.program_name
    prihdr['AUTHOR']  = ('Sebastian Schmidt')
    prihdr['DATE']    = (date_time.strftime("%A, %d. %B %Y %I:%M%p"), 'Date created')
    prihdr['ARRAY0']  = ('WAVE', 'Wavelength of composite spectra in Angstrom')
    for i in xrange(ncspec):
        # j counts the number of the array
        # 1, because ARRAY0 is set before
        j = 1 + i*4
        keyword1 = 'ARRAY' + str(j)
        keyword2 = 'ARRAY' + str(j+1)
        keyword3 = 'ARRAY' + str(j+2)
        keyword4 = 'ARRAY' + str(j+3)
        prihdr[keyword1] = ('FLUX', 'Flux of composite spectrum'+str(i))
        prihdr[keyword2] = ('EFLUX', 'Error of flux of composite spectrum'+str(i))
        prihdr[keyword3] = ('NORM', 'Normalized flux of composite spectrum'+str(i))
        prihdr[keyword4] = ('CONT', 'Continuum fit of composite spectrum'+str(i))

    prihdr['NCSPECS'] = (ncspec, 'Number of composite spectra')
    prihdr['Z_MAX']   = (settings.z_max, 'Maximal redshift binned')
    prihdr['Z_MIN']   = (settings.z_min, 'Minimal redshift binned')
    prihdr['Z_DELTA'] = (settings.z_delta, 'width of redshift bin')
    prihdr['NZBINS']  = (settings.z_bins, 'number of redshift bins')
#    prihdr['NSPECTRA']= (cspec.spectra_count, 'Total number of used / contributing QSOs')
#    prihdr['NSPEC']   = (nspec, 'Total number of QSOs considered')
    prihdr['VACUUM']  = (1, 'Wavelengths are in vacuum')
    prihdr['DC_FLAG'] = (1, 'Log-Linear Flag')
    prihdr['CRPIX1']  = (1, 'Starting pixel (1-indexed)')
    prihdr['COEFF0']  = np.log10(cspec.wave[0])
    prihdr['COEFF1']  = coeff1
    prihdr['CRVAL1']  = np.log10(cspec.wave[0])
    prihdr['CD1_1']   = coeff1
    prihdr['DR']      = (10, 'SDSS Data release used')

    print 'creation of arrays for HDU0'

    # zip arrays; will create one array of 3 tuples from the three
    # individual arrays
    print 'zipping arrays for HDU0'
    zipped_hdu = []
    zipped_hdu.append(compspectra[0].wave)
    for cspec in compspectra:
        zipped_hdu.append(cspec.flux)
        zipped_hdu.append(cspec.flux_error)
        zipped_hdu.append(cspec.normalized)
        zipped_hdu.append(cspec.powerlaw)

    print np.shape(zipped_hdu)
#    zipped_hdu = zip(zipped_hdu)
#    print np.shape(zipped_hdu)
    # currently wrong format, so we transpose the array
#    zipped_hdu = np.transpose(zipped_hdu)
    print np.shape(zipped_hdu)

    # create the Primary ouput HDU and print it to file
    hdu = fits.PrimaryHDU(zipped_hdu, header=prihdr)

    # Create 2nd HDU containing information on all QSOs
    # Calculate arrays, which will fill the table HDU
    # TODO: optimize creation of arrays! Lots of potential probably

    print 'creation of arrays for HDU1'
    counter_array        = np.zeros(ncspec)
    zem_index_array      = np.zeros(ncspec)
    zem_max_in_bin_array = np.zeros(ncspec)
    zem_min_in_bin_array = np.zeros(ncspec)

    for i, cspec in enumerate(compspectra):
        counter_array[i]        = cspec.spectra_count
        zem_index_array[i]      = cspec.zem_index
        zem_max_in_bin_array[i] = cspec.max_z_in_bin
        zem_min_in_bin_array[i] = cspec.min_z_in_bin
        if i % 500 == 0:
            print i, "spectra added to arrays"
    
    # write arrays to Table HDU:
    print 'write arrays to HDU1'
    TableHDU = fits.BinTableHDU.from_columns(
        fits.ColDefs([fits.Column(name='COUNT',  format='J', array=counter_array),
                      fits.Column(name='ZEM',    format='D', array=zem_index_array),
                      fits.Column(name='Z_MAX',  format='D', array=zem_max_in_bin_array),
                      fits.Column(name='Z_MIN',  format='D', array=zem_min_in_bin_array)]))

    # Now write new keywords to Table header
    table_header = TableHDU.header
    table_header['NCSPEC']  = (ncspec, 'Total number of composite spectra')

    # hdu.writeto('test.fits')
    print 'write FITS file'
    hdulist = fits.HDUList([hdu, TableHDU])
    # write to outfile and 'clobber=True' -> overwrite if existing
    hdulist.writeto(outfile, clobber=True)



