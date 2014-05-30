# Import used modules
# numpy: used for data storage etc.
import numpy as np
from pylab import *

# astropy is used to work with fits files
from astropy.io import fits
# re contains the function search, which is used to search for 
# spec and spSpec in the filenames
import re
# optimize contains the fitting function we use to fit the linear 
# function to the log-log spectrum data
from scipy import optimize
# obstools contains function to calculate E(B-V) based on Schlegel et. al dust maps
from astropysics import obstools
# contains astropy coordinate objects used to work with coordinates
from astropy.coordinates import ICRS, Galactic
from astropy import units as u
# contains function to call date and time
from datetime import datetime


# HDU numbers in this file are given starting from 0, since it is
# conform with astropy, although the FITS standard starts from 1!



######################################################################################
############################## Functions #############################################
######################################################################################


######################################################################################
############################## Functions working on files ############################
######################################################################################

# Need resolving power of spectrum?
def read_spec(qso, spec):
# Function whic is called to read spec FITS file

    #Open FITS file
    hdu = fits.open(qso)
    #Read in all necessary information from Header of
    spec.plateid = hdu[0].header['PLATEID']
    spec.fiberid = hdu[0].header['FIBERID']
    spec.MJD = hdu[0].header['MJD']

    # Coordinates:
    ra = hdu[0].header['PLUG_RA']
    dec = hdu[0].header['PLUG_DEC']
    spec.coordinates = ICRS(ra = ra, dec = dec, unit=(u.degree, u.degree))
    
    spec.beginwl = hdu[0].header['COEFF0']
    spec.deltawl = hdu[0].header['COEFF1']
    # Can't find cpix in header, assume it's 1.
    spec.cpix = 1
    spec.npix = hdu[1].header['NAXIS2']
    # Get the table data from HDU 1; it contains the flux
    # Change tabledata = hdu.data to tabledata = table.view(np.recarray)
    TableHDU1 = hdu[1].data
    spec.flux = TableHDU1['flux'].copy()
    spec.flux_error = TableHDU1['ivar'].copy()

    # Table data from HDU2; contains redshift, MJD and psf magnitudes
    TableHDU2 = hdu[2].data
    spec.z = float(TableHDU2['Z'])

    #magnitudes not needed for now. Did they cause memory overhead? 
    # TODO: check that!
#    spec.mag = TableHDU2['PSFMAG'].copy()

#TODO: Read in primary target information?


    # Call function which converts ra and dec to galactic coordinates

# assume logarithmic wavelength scale. Is used always I believe.?
    # run over all wavelengths (pixels) and set the wavelength array
    # as well as the status array and signal to noise ratio
    for i in xrange(spec.npix):
        spec.wave.append(10.0**(spec.beginwl + (i+1 - spec.cpix)*spec.deltawl))
        # check if flux_error > 0, because normalisation in snr divides by it
        if spec.flux_error[i] > 0: 
            spec.status.append(1)
            spec.snr.append(spec.flux[i] / spec.flux_error[i])
        # if flux_error > 0 simply set status and error to 0
        else:
            spec.status.append(0)
            spec.snr.append(0)
    
    del(TableHDU1)
    del(TableHDU2)
    hdu.close()


def read_spSpec(qso, spec):
# Function which is called to read spSpec FITS file

    #Open FITS file
    hdu = fits.open(qso)
    #Read in all necessary information from Header of HDU 0:
    spec.plateid = hdu[0].header['PLATEID']
    spec.fiberid = hdu[0].header['FIBERID']
    spec.MJD = hdu[0].header['MJD']
    spec.z = hdu[0].header['Z']
    spec.mag = hdu[0].header['MAG']

    # Coordinates:
    ra = hdu[0].header['RAOBJ']
    dec = hdu[0].header['DECOBJ']
    spec.coordinates = ICRS(ra = ra, dec = dec, unit=(u.degree, u.degree))

    spec.npix = hdu[0].header['NAXIS1']
    spec.beginwl = hdu[0].header['CRVAL1']
    spec.deltawl = hdu[0].header['CD1_1']
    spec.cpix = hdu[0].header['CRPIX1']
# Call a function which calculates galactic coordinates from ra and dec
    #Open ImageData from HDU 0 and thus retrieve the flux and the error
    ImageData = hdu[0].data.copy()
    spec.flux = ImageData[0, 0:spec.npix].copy()
    spec.flux_error = ImageData[1, 0:spec.npix].copy()
    
# assume logarithmic wavelength scale. Is used always I believe.?
    # run over all wavelengths (pixels) and set the wavelength array
    # as well as the status array and signal to noise ratio
    for i in xrange(spec.npix):
        spec.wave.append(10.0**(spec.beginwl + (i+1 - spec.cpix)*spec.deltawl))
        if spec.flux_error[i] > 0: 
            spec.status.append(1)
            spec.snr.append(spec.flux[i] / spec.flux_error[i])
        else:
            spec.status.append(0)
            spec.snr.append(0)

    del(ImageData)
    hdu.close()
# Do I need to divide ra by 15, like in C code? why is that done?
# What is primary target flag needed for?
    

def check_filetype(filename):
# Function which checks the filetype of a given file (string)
# simply searches for spSpec and spec in the filename 
    if re.search("spSpec", filename):
        return 1
    elif re.search("spec", filename):
        return 2
    else:
        return 0


################################################################################
############################## Dust corrections ################################
################################################################################

def get_Ebv(spec, dustmap):
# Function which calculates E(B-V) using obstools, based on Schlegel et. al 
# dustmaps.     
    # coordinates in galactic longitude, latitude
    l = spec.coordinates.galactic.l.deg
    b = spec.coordinates.galactic.b.deg
    
    spec.Ebv = obstools.get_SFD_dust(l, b, dustmap, interpolate=False)

def Gal_extinction_correction(spec):
# Function which applies the E(B-V) value to the corrections to the flux of
# the spectrum    
    # variables, taken from C program
    R_V = 3.08

    # A_V not necessary?
#    A_V = R_V*spec.Ebv

 #   A_lam = obstools.extinction_correction(spec.flux, spec.wave, spec.Ebv, Rv=R_V)
    print "flux: 1", spec.flux[123]
    spec.flux = obstools.extinction_correction(spec.flux, spec.wave, spec.Ebv, Rv=R_V)
#    print A_lam
#    for i in xrange(spec.npix):
        # check on A_lam? Maybe can return problematic value?
        # I suppose the next line is not necessary, as we already include Ebv and R_V into
        # the function?
        # A_lam *= A_V 
 
#        redfac = 10.0**(0.4*A_lam)
#        spec.flux[i] *= redfac
#        spec.flux_error[i] *= redfac
        
    print "flux: 2", spec.flux[123]
#    print "A_lam, redfac: ", A_lam, redfac


################################################################################
############################## Colors ##########################################
################################################################################

def create_colorcurves():
    # Open file into colorcurves and all lines into an array
    colorcurves = open("colorcurves" ,'r')
    lines = list(colorcurves)
    i = 0
    npix = 5763
    beginwl = 3.4742
    wavelength = []
    a = []
    for i in xrange(npix):
        wavelength.append(10.0**(beginwl + 0.0001*i))
    # Create arrays, which are to store the color curves
    ucurve = []
    gcurve = []
    rcurve = []
    icurve = []
    zcurve = []
    NU = NG = NR = NI = NZ = 0
    # Run over first few lines of the fiel in order to get the number
    # of elements, which will be stored in the curve arrays
    for num, line in enumerate(lines, 0):
        if "NU" in line:
            NU = eval(lines[num+1])
        if "NG" in line:
            NG = eval(lines[num+1]) 
        if "NR" in line:
            NR = eval(lines[num+1])
        if "NI" in line:
            NI = eval(lines[num+1])
        if "NZ" in line:
            NZ = eval(lines[num+1])
            break

    # Append an empty line at the end of the lines list, in order to
    # evaluate the last while query for the zcurve. Otherwise it would
    # be out of range, since i is set to i += 1.
    lines.append(" ")
    # run over all elements of the list, and enumerate the lines at
    # the same time
    for num, line in enumerate(lines, 0):
        # check if line contains the beginning of a color curve
        if "UCURVE" in line:
            # Set i to num+1, so that it is set to the line after UCURVE
            i = num+1
            # Now, read all following lines into curve array, until an
            # empty line follows
            while lines[i].strip() != '':
                # eval() is used to convert the strings found in each line
                # into floats
                ucurve.append(eval(lines[i]))
                i += 1
        # continue with file
        if "GCURVE" in line:
            i = num+1
            while lines[i].strip() != '':
                gcurve.append(eval(lines[i]))
                i += 1            
        if "RCURVE" in line:
            i = num+1
            while lines[i].strip() != '':
                rcurve.append(eval(lines[i]))
                i += 1
        if "ICURVE" in line:
            i = num+1
            while lines[i].strip() != '':
                icurve.append(eval(lines[i]))
                i += 1
        if "ZCURVE" in line:
            i = num+1
            while lines[i].strip() != '':
                zcurve.append(eval(lines[i]))
                i += 1

    # Now the curve arrays are defined and we need to calculate the actual color curves
    a = np.zeros((npix, 5))
    i = int(round((log10(ucurve[0][0]) - beginwl)/0.0001))
    # TODO: make the calculations look nicer by splitting into several lines
    for j in xrange(NU):
        while wavelength[i] < ucurve[j+1][0]:
            a[i][0] = (ucurve[j+1][1] - ucurve[j][1])/(ucurve[j+1][0] - ucurve[j][0])*(wavelength[i] - ucurve[j][0]) + ucurve[j][1]
            i += 1
    i = int(round((log10(gcurve[0][0]) - beginwl)/0.0001))
    for j in xrange(NG):
        while wavelength[i] < gcurve[j+1][0]:
            a[i][1] = (gcurve[j + 1][1] - gcurve[j][1]) / (gcurve[j+1][0]-gcurve[j][0]) * (wavelength[i] - gcurve[j][0]) + gcurve[j][1]
            i += 1
    i = int(round((log10(rcurve[0][0]) - beginwl)/0.0001))
    for j in xrange(NR):
        while wavelength[i] < rcurve[j+1][0]:
            a[i][2] = (rcurve[j + 1][1] - rcurve[j][1]) / (rcurve[j+1][0]-rcurve[j][0]) * (wavelength[i] - rcurve[j][0]) + rcurve[j][1]
            i += 1
    i = int(round((log10(icurve[0][0]) - beginwl)/0.0001))
    for j in xrange(NI):
        while wavelength[i] < icurve[j+1][0]:
            a[i][3] = (icurve[j + 1][1] - icurve[j][1]) / (icurve[j+1][0]-icurve[j][0]) * (wavelength[i] - icurve[j][0]) + icurve[j][1]
            i += 1
    i = int(round((log10(zcurve[0][0]) - beginwl)/0.0001))
    for j in xrange(NZ):
        while wavelength[i] < zcurve[j+1][0]:
            a[i][4] = (zcurve[j + 1][1] - zcurve[j][1]) / (zcurve[j+1][0]-zcurve[j][0]) * (wavelength[i] - zcurve[j][0]) + zcurve[j][1]
            i += 1

    #TODO: close colorcurves
    return a

def colors(spec, a):
# This function calculates the spectroscopic magnitudes based on the color curves
# calculated in create_colorcurves() and the psf magnitudes found in each
# spectrum's FITS file
    npix = 5763
    beginwl = 3.4742
    usum = gsum = rsum = isum = zsum = 0
    j = int(round((log10(spec.wave[0]) - beginwl)/0.0001))
    i = 0
    while i < spec.npix and j < npix:
        usum += a[j][0] * spec.flux[i]
        gsum += a[j][1] * spec.flux[i]
        rsum += a[j][2] * spec.flux[i]
        isum += a[j][3] * spec.flux[i]
        zsum += a[j][4] * spec.flux[i]
        i += 1
        j += 1
        
    if usum > 0 and gsum > 0 and rsum > 0 and isum > 0 and zsum > 0:
        spec.smag.append(22.5 - 2.5 * log10(usum))
        spec.smag.append(22.5 - 2.5 * log10(gsum))
        spec.smag.append(22.5 - 2.5 * log10(rsum))
        spec.smag.append(22.5 - 2.5 * log10(isum))
        spec.smag.append(22.5 - 2.5 * log10(zsum))
    else:
        spec.smag.append(0)
        spec.smag.append(0)
        spec.smag.append(0)
        spec.smag.append(0)
        spec.smag.append(0)

    # What to do with invalid values? 


################################################################################
############################## Powerlaw ########################################
################################################################################
    

def find_element_larger_in_arrays(wave, target_wave, npix):
    # Function which takes in a wavelength array wave and a target wavelength
    # It searches for the i-th element in wave, which is the first element
    # bigger than target_wave. 
    # Note: i is always larger than target_wave! If smaller than target_wave is 
    # wanted, take returnvalue - 1
    # TODO: Understand exactly how it works :) and think of a way to catch 
    # exception 
#    i = next(x[0] for x in enumerate(wave) if x[1] > target_wave)
    i = 0
    while True:
        if i < npix and wave[i] < target_wave:
            i += 1
        else:
            break
    return i

def fit_powerlaw(spec):
# This function fits the powerlaw to the spectrum, by taking the emission free regions
# emfree and fits a linear function to log-log flux- wavelength data
    emfree = []
    # define emission free regions
    emfree = np.array([[1280.0,1292.0],[1312.0,1328.0],[1345.0,1365.0],[1440.0,1475.0]])
    # we use 4 emission free regions
    emfree_regions_num = 4
    # create emtpy lists for the log wavelength, log flux and error data
    # TODO: think if lists should be numpy arrays instead!
    wave_log = []
    flux_log = []
    flux_error_log = []
    # set standard values for alpha, alpha_error, beta and delta
    spec.alpha = spec.alpha_error= spec.beta = spec.delta = -999
    # create a variable, which counts how many regions contain usable data
    # only fit a function if all 4 regions are usable
    emfree_regions_data = 0

    # Initialise spectrum's emfree matrix:
    spec.emfree = np.zeros((3,4))
    for i in xrange(emfree_regions_num):
        # calculate the array containing the emission free regions
        spec.emfree[0][i] = (1.0 + spec.z)*0.5*(emfree[i][0] + emfree[i][1])
        spec.emfree[1][i] = 0.0
        spec.emfree[2][i] = -1.0

        # Find the element in the wavelength array, which is larger / smaller 
        # than the beginning / end of the emission free region
        wave_em_interval_start = find_element_larger_in_arrays(spec.wave, (1.0 + spec.z)*emfree[i][0], spec.npix)
        # the -1 at the end takes into account, that the function always returns the bigger
        # value.
        # TODO: Check why siqr not the same as in C code and median neither. (slight abbreviations)
        wave_em_interval_end = find_element_larger_in_arrays(spec.wave, (1.0 + spec.z)*emfree[i][1], spec.npix) - 1
        # Calculate 75 and 25 percentile in order to calculate the semi-interquantile range
        percentile75 = np.percentile(spec.flux[wave_em_interval_start:wave_em_interval_end], 75)
        percentile25 = np.percentile(spec.flux[wave_em_interval_start:wave_em_interval_end], 25)
        siqr = (percentile75 - percentile25)/2
        median = np.percentile(spec.flux[wave_em_interval_start:wave_em_interval_end], 50)

        spec.emfree[1][i] = median
        spec.emfree[2][i] = siqr/(wave_em_interval_end - wave_em_interval_start)
        # if usable values, append emfree regions to log data arrays
        if spec.emfree[1][i] > 0 and spec.emfree[2][i] > 0:
            wave_log.append(log10(spec.emfree[0][i]))
            flux_log.append(log10(spec.emfree[1][i]))
            flux_error_log.append(log10(1.0 + spec.emfree[2][i]/spec.emfree[1][i]))
            # count region als usable
            emfree_regions_data += 1
        
    # Fit a linear function to the log data:
    # define our (line) fitting function
    def func(x, a, b):
        return a*x + b
    # only continue, if 4 regions contain usable data
    if emfree_regions_data == 4:
        # fit linear function func to out log arrays
        # coeff: fitting parameters
        # pcov: covariance matrix, used to retrieve error for alpha.
        coeff, pcov = optimize.curve_fit(func, wave_log, flux_log, p0=(-2,5), sigma=flux_error_log)

        # assign coefficients to our spectrum
        spec.beta = coeff[0]
        spec.alpha = -spec.beta - 2
        spec.delta = coeff[1]
        try:
            spec.alpha_error = float(sqrt(pcov[0][0]))
        except TypeError:
            print "Fitting problem"

        # use the fitted coefficients to calculate the powerlaw continuum
        for i in xrange(spec.npix):
            spec.powerlaw.append(10.0**(coeff[1] + coeff[0]*log10(spec.wave[i])))

    # if we don't have 4 usable regions, set everything to 0
    else:
        spec.beta = 0
        spec.alpha = 0
        spec.delta = 0
        spec.alpha_error = 0
        for i in xrange(spec.npix):
            spec.powerlaw.append(0)    
            # Think about if 0 is a good value (currently checked in build_compspec)

    
    del(spec.emfree)
    # TODO: values differ slightly from c program!
    # probably because different fitting algorithm
    
################################################################################
############################## Compspec ########################################
################################################################################
    
def build_compspec(cspec, spec):
# this function calculates the 'composite spectrum' from the individual spectra
# by adding flux/powerlaw values in a certain Restframe wavelength range
    # Define Restrange array:
    Restrange = []
    Restrange = np.array([[1095,1150],[982,1010],[955,968],[912,947]])

    # Set the interval bounds of the QSO wavelengths, which are used
    if (spec.z + 1)*Restrange[0][0] > spec.wave[0]:
        spec.lambmin = (spec.z + 1)*Restrange[0][0]
    else:
        spec.lambmin = spec.wave[0]
    spec.lambmax = (spec.z + 1)*Restrange[0][1]

    # if first wavelength of spectrum bigger than lambmax, discard QSO
    if spec.lambmax < spec.wave[0]:
        return 9

    # First: align cspec arrays with spec arrays:
    # iterator and citerator are the number of elements which seperate the arrays
    # of compspec and the spectrum. They are added to the actual iterator, of
    # the for loop, in order to align the wavelength, ... arrays, so that we
    # actually add the correct contributions of each spectrum to the correct 
    # wavelength in the composite spectrum    
    iterator = find_element_larger_in_arrays(spec.wave, cspec.wave[0], spec.npix)
    citerator = find_element_larger_in_arrays(cspec.wave, spec.wave[0], spec.npix)
    # Run over all pixels / wavelengths of the spectrum
    # - iterator, because we don't want to access elements outside of array bounds
    # from 1 to pixels - iterator - 1, because we check for lambmin and lambmax
    # by looking at element i - 1, i + 1 respectively
    for i in xrange(1, spec.npix - iterator - 1):
        # if a lot of checks are true, 
        if(spec.wave[i+iterator]                             != 0             and
           spec.wave[i+iterator] + spec.wave[i+iterator - 1] > 2*spec.lambmin and
           spec.wave[i+iterator] + spec.wave[i+iterator + 1] < 2*spec.lambmax and
           isnan(spec.flux[i+iterator])                      == 0             and
           isnan(spec.flux_error[i+iterator])                == 0             and
           spec.flux_error[i+iterator]                       != 0             and
           spec.powerlaw[i+iterator]                         != 0):
            # we actually use the spectrum to add to the composite spectrum
            cspec.sum[i+citerator] += spec.flux[i+iterator]/spec.powerlaw[i+iterator]
            cspec.sum2[i+citerator] += (spec.flux[i+iterator]/spec.powerlaw[i+iterator])**2
            # and we count this contribution to this wavelength
            cspec.nhist[i+citerator] += 1
        else: 
            # if checks not met, go to next pixel
            continue
        # if we would be outside of lambmax in the next iteration, stop for loop
        if spec.wave[i+iterator] + spec.wave[i+iterator + 1] > 2*spec.lambmax: 
            break
            
    # if everything goes well and we used the QSO, return 0
    return 0

def combine_compspecs(compspectra):
# A function which combines all compsectra from each process into one
    num_compspec = len(compspectra)
    npix = 5763
    compspec_sum = comp_spectrum(npix)
    # Run over all compspectra and simply add the values for each pixel
    # of sum, sum2 and nhist as well as calculate the total number of
    # contributing QSOs from the spectra_counts
    for compspec in compspectra:
        for j in xrange(npix):
            compspec_sum.sum[j] += compspec.sum[j]
            compspec_sum.sum2[j] += compspec.sum2[j]
            compspec_sum.nhist[j] += compspec.nhist[j]
        compspec_sum.spectra_count += compspec.spectra_count
    
    return compspec_sum


def statistics(cspec, spec):
# a function which calculates certain values
    # hardcoded the pixel range of 3000, same as in C program. 
    # TODO: check if can be done nicer
    # run over pixel range
    for i in xrange(3000):
        # if sum, nhist and sum2 are good, 
        if cspec.sum[i] > 0 and cspec.nhist[i] > 0 and cspec.sum2[i] > 0:
            # we add the flux of the composite spectrum as:
            # sum / nhist
            cspec.flux.append(cspec.sum[i]/cspec.nhist[i])
            cspec.flux_error.append(cspec.sum2[i]/cspec.nhist[i] - cspec.flux[i]**2)
            if cspec.nhist[i] > 1:
                cspec.flux_error[i] = sqrt(cspec.flux_error[i]/(cspec.nhist[i] - 1))
            else:
                cspec.flux_error[i] = 0
        else:
            # if bad, we just append 0.
            cspec.flux.append(0)
            cspec.flux_error.append(0)

    # calculate mean_a, sigma_a, median_a and siqr_a
    nspec = len(spec)
#    print nspec, spec
    # Create an array containing all alpha values of spectra, to which a
    # powerlaw could be fitted. Others don't contribute.
    alpha_array = np.zeros(nspec)
    for i in xrange(nspec):
        if spec[i].alpha > -900:
            alpha_array[i]  = spec[i].alpha

    # use alpha_array to calculate quantities
    percentile75   = np.percentile(alpha_array, 75)
    percentile25   = np.percentile(alpha_array, 25)
    cspec.siqr_a   = (percentile75 - percentile25)/2
    cspec.median_a = np.percentile(alpha_array, 50)
    cspec.mean_a   = np.mean(alpha_array)
    cspec.sigma_a  = np.std(alpha_array)

  

################################################################################
############################## Output ##########################################
################################################################################

def build_fits_file(cspec, spec, outfile, settings):
# Build the actual FITS output file
    # set a few variables
    i = 0
    coeff1 = 0.0001
    while cspec.nhist == 0:
        i += 1
    sidx = i
    npix = 0
    date_time = datetime.now()

    hdu0_row1 = []
    hdu0_row2 = []
    hdu0_row3 = []
    for i in xrange(sidx, 3000):
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
    prihdr['NSPECTRA']= (cspec.spectra_count, 'Total number of contributing QSOs')
    prihdr['NSPEC']   = (len(spec), 'Total number of contributing QSOs')
    prihdr['VACUUM']  = (1, 'Wavelengths are in vacuum')
    prihdr['DC_FLAG'] = (1, 'Log-Linear Flag')
    prihdr['CRPIX1']  = (1, 'Starting pixel (1-indexed)')
    prihdr['COEFF0']  = log10(cspec.wave[sidx])
    prihdr['COEFF1']  = coeff1
    prihdr['CRVAL1']  = log10(cspec.wave[sidx])
    prihdr['CD1_1']   = coeff1
    prihdr['DR']      = (10, 'SDSS Data release used')


    # run over all pixels and append elements to new arrays
    for i in xrange(npix):
        hdu0_row1.append(1 - cspec.flux[sidx+i])
        hdu0_row2.append(cspec.flux_error[sidx+i])
        hdu0_row3.append(cspec.nhist[sidx+i])
    # zip arrays; will create one array of 3 tuples from the three
    # individual arrays
    zipped_hdu0 = zip(hdu0_row1, hdu0_row2, hdu0_row3)
    # currently wrong format, so we transpose the array
    zipped_hdu0 = np.transpose(zipped_hdu0)

    # create the Primary ouput HDU and print it to file
    hdu = fits.PrimaryHDU(zipped_hdu0, header=prihdr)

    # Create 2nd HDU containing information on all QSOs
    # Calculate arrays, which will fill the table HDU
    mjd_array         = map(lambda spec: spec.MJD, spec)
    plate_array       = map(lambda spec: spec.plateid, spec)
    fiber_array       = map(lambda spec: spec.fiberid, spec) 
    alpha_array       = map(lambda spec: spec.alpha, spec)
    alpha_error_array = map(lambda spec: spec.alpha_error, spec)
    ra_array          = map(lambda spec: spec.coordinates.ra.deg, spec)
    dec_array         = map(lambda spec: spec.coordinates.dec.deg, spec)
    l_array           = map(lambda spec: spec.coordinates.galactic.l.deg, spec)
    b_array           = map(lambda spec: spec.coordinates.galactic.b.deg, spec)
    zem_array         = map(lambda spec: spec.z, spec)
    smag_array        = map(lambda spec: spec.smag, spec)
    #flag_array 

    # write arrays to Table HDU:
    TableHDU = fits.new_table(
        fits.ColDefs([fits.Column(name='MJD',    format='D', array=mjd_array),
                      fits.Column(name='PLATE',  format='D', array=plate_array),
                      fits.Column(name='FIBER',  format='D', array=fiber_array),
                      fits.Column(name='RA',     format='D', array=ra_array),
                      fits.Column(name='DEC',    format='D', array=dec_array),
                      fits.Column(name='ALPHA',  format='D', array=alpha_array),
                      fits.Column(name='EALPHA', format='D', array=alpha_error_array),
                      fits.Column(name='l',      format='D', array=l_array),
                      fits.Column(name='b',      format='D', array=b_array),
                      fits.Column(name='ZEM',    format='D', array=zem_array),
                      fits.Column(name='MAGSPEC',format='PD()', array=smag_array)]))
     #                 fits.Column(name='FLAG',   format='E', array=flag_array)]))
    # TODO: check if I can change smag column a little bit


    # Now write new keywords to Table header
    table_header = TableHDU.header
    table_header['NSPEC']  = (len(spec), 'Total number of Quasars')
    table_header['MEAN_A'] = (cspec.mean_a, 'Mean spectral index')
    table_header['SIGMA_A']= (cspec.sigma_a, 'Standard deviation on alpha')
    table_header['MED_A']  = (cspec.median_a, 'Median spectral index')
    table_header['SIQR_A'] = (cspec.siqr_a, '68% semi-interquartile range on alpha')

    

    # hdu.writeto('test.fits')
    hdulist = fits.HDUList([hdu, TableHDU])
    # write to outfile and 'clobber=True' -> overwrite if existing
    hdulist.writeto(outfile, clobber=True)



################################################################################
############################## Help screen #####################################
################################################################################


def help():
    print "####################    SDSScompspec in Python    ####################"
    print ""
    print "Usage:"
    print "./PyS_SDSScompspec <input list file> --options"
    print "possible options:"
    print "    --dust:    activates galactic dust corrections"
    print "    -h, --help: print this help"


################################################################################
############################## Classes #########################################
################################################################################
 
# Class for the spectra
class spectrum:
    def __init__(self):
        self.beginwl     = 0
        self.deltawl     = 0
        self.lambmin     = 0
        self.lambmax     = 0
        self.cpix        = 0
        self.z           = 0
        # Fitting parameters
        self.alpha       = 0
        self.alpha_error = 0
        self.beta        = 0
        self.delta       = 0
        # Coordinates
        # using astropy.coordinates object. Currently using v0.3.2 of astropy
        # from v0.4 (in develeopment atm) on exchanged with SkyCoord. 
        self.coordinates = ICRS(ra = 0, dec = 0, unit=(u.degree, u.degree))

        self.npix           = 0
        self.flux           = []
        self.flux_error     = []
        self.powerlaw       = [] 
        self.powerlaw_error = []
        self.wave           = []
        self.emfree         = []         # Emission free regions array
        self.mag            = []
        self.mag_error      = []
        self.smag           = []
        self.status         = []
        self.snr            = []

        self.filename = []
        self.plateid = 0
        self.fiberid = 0
        self.MJD = 0

# Class for the composite spectrum
class comp_spectrum:
    def __init__(self, npix):
        self.wave          = []
        self.flux          = []
        self.flux_error    = []
        self.sum           = []
        self.sum2          = []
        self.nhist         = []
        self.spectra_count = 0
        self.mean_a        = 0
        self.sigma_a       = 0
        self.median_a      = 0
        self.siqr_a        = 0
        for i in xrange(npix):
            self.sum.append(0)
            self.sum2.append(0)
            self.nhist.append(0)
            self.wave.append(10.0**(3.58020 + 0.0001*i))

# Class for the program settings
class program_settings:
    def __init__(self):
        self.dust         = 0
        self.nprocs       = 0
        self.program_name = 'PyS_SDSScompspec'
        
