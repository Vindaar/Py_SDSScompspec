# This module contains different functions to perform flux correction as
# well as dust correction

import numpy as np

# fitting in flux correction
from scipy.integrate import quad

# obstools contains function to calculate E(B-V) based on Schlegel et. al dust maps
from astropysics import obstools

# contains astropy coordinate objects used to work with coordinates
from astropy.coordinates import SkyCoord
from astropy import units as u

# fitsio based on C
import fitsio

from SDSSmodules.SDSSfiles import *

# cython compiled dust_extinction routine
# Dust extinction routine calculates the actual values for the dust correction. 
# check whether dust_extinction_array was compiled in cython or not. 
# If not, compile now
# TODO: check real speed difference between python and cython
# factor ~60 faster
try:
    from dust_extinction_array import dust_extinction_array
except ImportError:
    print 'It seems like dust_extinction_array.pyx was not compiled properly.'
    print 'This is done by typing:'
    print 'python setup.py build_ext --inplace'
    print 'in the folder of PyS_SDSScompspec.py.'
    print ''
    if raw_input("Do you wish to compile it now? (y/N) ") == ('y' or 'Y'):
        import os
        if os.system("python setup.py build_ext --inplace") == 0:
            try:
                from dust_extinction_array import dust_extinction_array
                print 'compilation successful!'
                print 'starting program...'
            except ImportError:
                print 'compilation not successful'
                import sys
                sys.exit("Exit program")
        else:
            import sys
            sys.exit("setup.py is missing or some unknown error occured")
    else:
        print "Do you wish to use the python based dust_extinction_array instead?"
        print "WARNING: it is considerably slower!"
        if raw_input("(y/N) ") == ('y' or 'Y'):
            from dust_extinction_python import dust_extinction_array_python as dust_extinction_array 
        else:
            import sys
            sys.exit("Exit program")


#####################################################################################
################################### Flux correction #################################
#####################################################################################


def partial_water_pressure(P, T, h):
    # calculate partial pressure of water
    # 1. calculate equilibrium vapor pressure of water at pressure P and temp T
    # based on Buck 1981
    # http://www.public.iastate.edu/~bkh/teaching/505/arden_buck_sat.pdf 24.09.14
    # 2. relative Hum defined as: h = ppW / eppW * 100
    # where ppW is the partial pressure of water in the gas, eppW the equilibrium vapor 
    # pressure of water at that pressure and temp
    if np.isnan(P) or np.isnan(T) or np.isnan(h):
        print P, T, h
        import sys
        sys.exit('ERROR in partial_water_pressure(): Pressure, temperature and humidity values all NAN. Something bad happened before!!')
    eppW = (1.0007 + 3.46 * 10**(-6) * P) * (6.1121)*np.exp(17.502*T / (240.97 + T))
    ppW  = h * eppW / 100
    #print 'ppW:', ppW, eppW, P
    return ppW


def n(x, P, ppW, T):
    # function to caluclate refractive index of air
    # refractivity of atmosphere based on wavelength, pressure, partial water pressure and temperature
    # Marini, J.W., NASA Technical Report X-591-73-351, (1973). 
    # http://physics.ucsd.edu/~tmurphy/apollo/doc/MM.pdf 24.09.14
    #print 'propert', x, P, ppW, T
    refr = (287.604 + 1.6288/(x**2) + 0.0136/(x**4)) * (P / 1013.25) * (1 / (1 + 0.003661*T)) - 0.055 * (760 / 1013.25) * (ppW / (1 + 0.00366 * T))
    # calculate refractive index n from refractivity refr
    # refr = 10^6 (n - 1)
    val = refr * 10**(-6) + 1
    return val

def calc_delta_y(lam, reference_lam, P, ppW, T, Z):
    # function which calculates the shift in arcseconds of the flux maximum
    # based on the refractive index of the atmosphere relative to a given
    # wavelength
    val = 206265 * ( n(lam, P, ppW, T) - n(reference_lam, P, ppW, T) ) * np.tan(Z)
    return val


def perform_flux_corr_ind_exp(spec, settings, resid_corr):

    # First start with residual correction for Balmer problem
    # resid_corr array starts from wavelength 10**(3.5496)
    if np.size(resid_corr) > 1:
        start_resid_corr = int((spec.beginwl - 3.5496)*1000)
        resid_corr = resid_corr[start_resid_corr:start_resid_corr+spec.npix]
    spec.flux       = spec.flux       / resid_corr
    spec.flux_error = spec.flux_error / resid_corr

    # Now perform correction for loss of flux due to offset in position of fiber

    # fiber diameter 2 arcseconds on sky (120 microns)
    # fiber radius in arcseconds
    a = 1
    # lam = wavelength in microns
    # Angstrom : 10e-10 m, micron : 10e-6 m
    lam = spec.wave * 10**(-10) / (10**(-6))
    reference_lam = 4000*10**(-4)
    print 'ghe'

    seeing_arr, altitude_arr, P_arr, T_arr, humidity_arr, snr_arr, flux_list, error_list, wave_list, loglam = get_array_from_ind_exposures(spec, settings, return_array = 1)
    ppW_arr = partial_water_pressure(P_arr, T_arr, humidity_arr)

    def flux_correction(delta_y):
        #print delta_y
        func = lambda x: 1/(2*np.pi)*(1 - np.exp(-(a**2 - delta_y**2 + 2*delta_y*np.sin(x)*(delta_y*sin(x) + np.sqrt(a**2 + delta_y**2*(np.sin(x)**2 - 1))))/(2*sigma**2)))
        flux_corr = quad(func, 0, 2*np.pi)
        return flux_corr[0]

    def flux_base():
        #print delta_y
        func = lambda x: 1/(2*np.pi)*(1 - np.exp(-(x)/(2*sigma**2)))
        flux_corr = quad(func, 0, 2*np.pi)
        return flux_corr[0]        
        

    flux_corr = []
    num_exp = np.size(seeing_arr)
    for i in xrange(num_exp):
        # zenith angle in radian
        P        = P_arr[i]
        # Convert pressure from Inches of Hg to HPa
        P = P * 3386.389 / 100
        T        = T_arr[i]
        ppW      = ppW_arr[i]
        seeing   = seeing_arr[i]
        altitude = altitude_arr[i]
        snr      = snr_arr[i]
        flux     = flux_list[i]
        error    = error_list[i]
        lam      = 10**(wave_list[i])*10**(-4)

        Z = (90 - altitude)*2*np.pi / 360
        # sigma from seeing
        sigma = seeing / (2*np.sqrt(2 * np.log(2)))

        delta_y = calc_delta_y(lam, reference_lam, P, ppW, T, Z)
        print lam[0:10], P, ppW, T, Z
        print delta_y[0:25]

        # hopefully this works!

        #flux_5400 = flux_correction(calc_delta_y(4000*10**(-4), 5400*10**(-4), P, ppW, T))
        #print 'flux_base', flux_5400
        flux_base = flux_correction(0)
        flux_qso  = empty(np.size(delta_y))
        print np.size(delta_y)
        for i in xrange(np.size(delta_y)):
            flux_qso[i]  = flux_correction(delta_y[i])
            #print delta_y[i], i, flux_qso[i]
        #flux_qso = flux_correction(delta_y)
        flux_corr.append(flux_qso / flux_base)
        print np.shape(flux_corr), np.shape(spec.flux)

    #flux_corr_blue = np.empty(num_exp/2)
    #flux_corr_red  = np.empty(num_exp/2)

    flux_corr_final = np.empty(np.size(loglam))
    flux_corr_final.fill(1.0)

    for i in xrange(np.size(loglam)):
        wave_pixel = loglam[i]
        num = 0
        flux = []
        snr  = []
        for j, list in enumerate(wave_list):
            #print list[i]
            contained = np.isclose(list, wave_pixel, atol=1e-5)
            if np.any(contained):
                contained = np.where(contained)[0][0]
                print j, contained
                flux.append(flux_corr[j][contained])
                snr.append(snr_arr[j])
        print 'shape', np.shape(flux), np.shape(snr)
        flux_corr_final[i] = np.ma.average(flux, weights=snr)
        #num = np.size(np.where(np.isclose(wave_list[i][:], wave_pixel, atol=1e-5))[0])
        print 'wave_pixel', wave_pixel
        #print np.around(wave_list[0], decimals=4)
        #print wave_list[0]
        print num


    # wave_list is in loglam!


    # for i in xrange(num_exp/2 - 1):
    #     if i == 0:
    #         index = np.in1d(wave_list[0], wave_list[1])
    #         #wave_common = np.intersect1d(wave_list[0], wave_list[1])
    #     else:
    #         index = np.in1d(wave_list[0][index], wave_list[i+1])
    #         #wave_common = np.intersect1d(wave_common, wave_list[i+1])

    # wave_common = wave_list[0][index]
    # print 'number common elements', np.size(index)
    # print np.size(np.unique(wave_list[0]))

    # print 'shape:', np.shape(flux_corr), num_exp
    # flux_corr_blue = np.ma.average(flux_corr[:num_exp/2-1], axis=1, weights=snr_arr[:num_exp/2-1])
    # flux_corr_red  = np.ma.average(flux_corr[num_exp/2:], axis=1, weights=snr_arr[num_exp/2:])

    # start_overlap  = find_element_larger_in_array(wave_list[0], wave_list[num_exp/2][0], np.size(wave_list[0]))
    # end_overlap    = find_element_larger_in_array(wave_list[num_exp/2], wave_list[0][-1], np.size(wave_list[0]))
    # flux_corr_final = flux_corr_blue[0:start_overlap-1]
    # overlap = (flux_corr_blue[start_overlap:], flux_corr_red[:end_overlap])
    # overlap_mean = np.ma.average(overlap, axis=1)
    # flux_corr_final.append(overlap_mean)
    # flux_corr_final.append(flux_corr_red[end_overlap:])
    # print np.size(flux_corr_final)
        
    return flux_corr_final

    #flux_correction = dblquad(lambda r, theta: (1 / sigma * np.sqrt(2*np.pi))*np.exp( -r**2 / 2*sigma**2), 0, 2*np.pi, lambda x: 0, lambda x: np.sqrt(a**2 - delta_y**2 * np.sin(theta)) - delta_y * np.cos(theta))



def perform_flux_correction(spec, settings, resid_corr, seeing='seeing50', check_weather = 0):
    # First start with residual correction for Balmer problem
    # resid_corr array starts from wavelength 10**(3.5496)
    if np.size(resid_corr) > 1:
        start_resid_corr = int((spec.beginwl - 3.5496)*1000)
        resid_corr_temp = resid_corr[start_resid_corr:start_resid_corr+spec.npix]
    spec.flux       = spec.flux       / resid_corr_temp
    spec.flux_error = spec.flux_error / resid_corr_temp

    # Now perform correction for loss of flux due to offset in position of fiber

    # fiber diameter 2 arcseconds on sky (120 microns)
    # fiber radius in arcseconds
    a = 1
    # lam = wavelength in microns
    # Angstrom : 10e-10 m, micron : 10e-6 m
    lam = spec.wave * 10**(-10) / (10**(-6))
    #reference_lam = 5400*10**(-4)
    from SDSSmodules.SDSSfiles import get_array_from_ind_exposures

    #try:
    seeing50, seeing80, altitude, P, T, ppW = get_array_from_ind_exposures(spec, settings, check_weather=check_weather)
    #except ValueError:
    #    raise ValueError('No weather data was found!')
    spec.altitude = altitude
    def n(x, P, ppW, T):
        # refractivity of atmosphere based on wavelength, pressure, partial water pressure and temperature
        # Marini, J.W., NASA Technical Report X-591-73-351, (1973). 
        # http://physics.ucsd.edu/~tmurphy/apollo/doc/MM.pdf 24.09.14
        #print 'propert', x, P, ppW, T
        refr = (287.604 + 1.6288/(x**2) + 0.0136/(x**4)) * (P / 1013.25) * (1 / (1 + 0.003661*T)) - 0.055 * (760 / 1013.25) * (ppW / (1 + 0.00366 * T))
        # calculate refractive index n from refractivity refr
        # refr = 10^6 (n - 1)
        val = refr * 10**(-6) + 1
        return val

    # zenith angle in radian
    Z = (90.0 - altitude)*2*np.pi / 360.0
    # sigma from seeing
    #seeing = 0.1
    if seeing == 'seeing50':
        sigma = seeing50 / (2*np.sqrt(2 * np.log(2)))
    elif seeing == 'seeing80':
        sigma = seeing80 / (2*np.sqrt(2 * np.log(2)))
    #print 2*np.sqrt(2*np.log(2))
    #print 'seeing is: ', seeing, sigma
    #sigma = sigma * 
    
    # test sigma influence
    #sigma = 0.2

    def calc_delta_y(lam, reference_lam, P, ppW, T, Z):
        val = 206265 * ( n(lam, P, ppW, T) - n(reference_lam, P, ppW, T) ) * np.tan(Z)
        return val

    # delta_y = calc_delta_y(lam, reference_lam, P, ppW, T, Z)
    # hopefully this works!
    def flux_correction(delta_y):
        #print delta_y
        #R = np.sqrt((a**2 - delta_y**2 + 2*delta_y*np.sin(0)*(delta_y*np.sin(0) + np.sqrt(a**2 + delta_y**2*(np.sin(0)**2 - 1)))))
        #print 'R', R, sigma, a, delta_y
        func = lambda x: 1/(2*np.pi)*(1 - np.exp(-(a**2 - delta_y**2 + 2*delta_y*np.sin(x)*(delta_y*np.sin(x) + np.sqrt(a**2 + delta_y**2*(np.sin(x)**2 - 1))))/(2*sigma**2)))
        flux_corr = quad(func, 0, 2*np.pi)
        return flux_corr[0]

    def flux_base():
        #print delta_y
        # is weird
        func = lambda x: 1/(2*np.pi)*(1 - np.exp(-(a**2)/(2*sigma**2)))
        flux_corr = quad(func, 0, 2*np.pi)
        return flux_corr[0]        


    #flux_base1 = flux_correction(calc_delta_y(4000*10**(-4), 5400*10**(-4), P, ppW, T, Z))
    #print 'flux_base', flux_5400
    reference_lam = 5400*10**(-4)
    delta_y = calc_delta_y(lam, reference_lam, P, ppW, T, Z)
    delta_y_5400 = delta_y
    print 'delta y should be :', delta_y

    # change radius of fiber artificially to check calc
    #a = 2*sigma
    #print 'a', a

    flux_base1 = flux_correction(0)
    flux_base2 = flux_base()
    #print 'delta 0', flux_base1, flux_base2

    print 'delta', min(delta_y), max(delta_y)
    flux_5400  = np.empty(np.size(delta_y))
    for i in xrange(np.size(delta_y)):
        flux_5400[i]  = flux_correction(delta_y[i])
        # if flux_5400[i] < 0.9:
        # #if flux_5400[i] < 0.999 and flux_5400[i] > 0.9985:
        #    print 'flux5400', flux_5400[i], delta_y[i]


    reference_lam = 4000*10**(-4)
    delta_y = calc_delta_y(lam, reference_lam, P, ppW, T, Z)
    delta_y_4000 = delta_y
    flux_4000  = np.empty(np.size(delta_y))
    # the actual integratoion is done for all delta_y individually. That means we do
    # about 4500 integrations for each spectrum
    # in total 166500 * 4500 ~ 750 mio integrations. Code is very slow...!
    # write in C?
    for i in xrange(np.size(delta_y)):
        flux_4000[i]  = flux_correction(delta_y[i])
        #if flux_qso[i] < 0.9:
        #if flux_qso[i] < 0.999 and flux_qso[i] > 0.9985:
            #print 'flux4000', flux_qso[i], delta_y[i]
        #print delta_y[i], i, flux_qso[i]
    #flux_qso = flux_correction(delta_y)
    #flux_corr = flux_qso * flux_5400) #/ flux_base1# / flux_base2
    #flux_qso / flux_base2
    #flux_5400 = flux_5400 / flux_base1
    #flux_qso = flux_qso / flux_base1

    # with the returned arrays, what we're going to do is the following:
    # multiply the spectrum by flux_5400 to revert the 'wrongly calibrated spectrum' to 
    # something resembling an uncalibrated spectrum
    # and then divide by flux_4000 to get a properly corrected spectrum
    return flux_4000, flux_5400, flux_base2


    #flux_correction = dblquad(lambda r, theta: (1 / sigma * np.sqrt(2*np.pi))*np.exp( -r**2 / 2*sigma**2), 0, 2*np.pi, lambda x: 0, lambda x: np.sqrt(a**2 - delta_y**2 * np.sin(theta)) - delta_y * np.cos(theta))


def apply_flux_correction(spec, correction_fits_file, settings, resid_corr):

    try:
        hdu = fitsio.FITS(correction_fits_file.rstrip())
    except ValueError:
        return 1
    #print correction_fits_file.rstrip()
    flux_corr = hdu[1]['corr'][:]
    spec.flux = spec.flux / flux_corr
    #print flux_corr
    hdu.close()
    return 0


def apply_resid_corr(spec, resid_corr):
    # Function not needed, if we apply the full flux correction. But for analysis purposes
    # it might be interesting to only apply the resid correction, to see the effect it has
    if np.size(resid_corr) > 1:
        start_resid_corr = int((spec.beginwl - 3.5496)*1000)
        resid_corr_temp = resid_corr[start_resid_corr:start_resid_corr+spec.npix]
    spec.flux       = spec.flux       / resid_corr_temp
    spec.flux_error = spec.flux_error / resid_corr_temp


################################################################################
############################## Dust corrections ################################
################################################################################

def fill_coordinate_arrays_from_buffer(coords, spec, dustmap, i, j):
# a function which fills the arrays of the coordinate_arrays
# object from the buffer, i.e. from index i to index i+j
# after that get Ebv values from dustmap
    coords.ra_array[i:i+j]  = map(lambda spec: spec.coordinates.ra.deg,  spec[i:i+j])
    coords.dec_array[i:i+j] = map(lambda spec: spec.coordinates.dec.deg, spec[i:i+j])
    temp_sky_objs    = SkyCoord(ra = coords.ra_array[i:i+j]*u.degree, dec = coords.dec_array[i:i+j]*u.degree, frame='fk5')
    coords.l_array[i:i+j]   = temp_sky_objs.galactic.l.deg
    coords.b_array[i:i+j]   = temp_sky_objs.galactic.b.deg
    Ebv = obstools.get_SFD_dust(coords.l_array[i:i+j], coords.b_array[i:i+j], dustmap, interpolate=2)
    k = 0
    for i in xrange(i, i+j):
        spec[i].Ebv = Ebv[k]
        k += 1

def Gal_extinction_correction(spec):
# Function which applies the E(B-V) value to the corrections to the flux of
# the spectrum    
    # variables, taken from C program
    R_V = 3.08

    # A_V not necessary?
    A_V = R_V*spec.Ebv

#    for i in xrange(spec.npix):
#        A_lam = obstools.extinction_correction(None, spec.wave[i], spec.Ebv, Rv=R_V)
#        A_lam = dust_extinction(wave[i], R_V)
#        A_lam *= A_V
#        redfac = 10**(0.4*A_lam)
#        flux[i]       *= redfac
#        flux_error[i] *= redfac

# version using obstools
#    spec.flux = obstools.extinction_correction(spec.flux, spec.wave, spec.Ebv, Rv=R_V)
    spec.flux, spec.flux_error = dust_extinction_array(spec.flux, spec.flux_error, spec.wave, R_V, A_V, spec.npix)






    

