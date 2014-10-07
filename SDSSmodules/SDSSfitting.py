# This module contains all the differnt continuum fitting functions

import numpy as np

# optimize contains the fitting function we use to fit the linear 
# function to the log-log spectrum data
from scipy import optimize
from linfit import linfit

# copy object
from copy import deepcopy

from SDSSutilities import *


################################################################################
############################## Fitting Powerlaw + CompSpec #####################
################################################################################



def fit_powerlaw_times_basespec(spec, settings, compspec_flux = 0, return_data = 0, deviation_factor = 3.0, zem_index = -999):
# This function fits the powerlaw to the spectrum, by taking the emission free regions
# emfree and fits a linear function to log-log flux- wavelength data
# In contrast to fit_powerlaw() it doesn't take the median of the intervals, but actually
# works on the individual data points in the intervals

# TODO: properly comment function, important since some things not so obvious
    # define emission free regions
#    emfree = np.array([[1280.0,1292.0],[1312.0,1328.0],[1345.0,1365.0],[1440.0,1475.0], [1610,1790]])

    zem     = spec.z
    z_max   = settings.z_max
    z_min   = settings.z_min
    z_delta = settings.z_delta
#    emfree = np.array([[1445.0, 1455.0], [1973.0, 1983.0]])
    
    # currently uselessly two regions, because code doesn't work with one
    # however, don't want to restructure, because eventually probably not the smartest
    # thing to take the whole region!
    emfree = np.array([[1205.0, 1400.0], [1400.0, 1800.0]])

    #    hdu = fitsio.FITS('CalcCompFits/DR10_noBALetc_newVarFits.fits')
    # hdu = fitsio.FITS('./CalcCompFits/dr9_compspectrum_2.5-3.fits')

    # calculate base spectrun
    # seperate smoothing of ly alpha forest and rest of spectrum

    # idea:
    # to filter ly alpha forest. multiply flux by s/n of pixel
    # -> pixel with low values no effect. 
    # if S/N in general low, effect small, 
    # at end, divide by average s/n to get into same magnitude again

    #print nsize
    flux_means = smooth_array(spec.flux, 3, 1/spec.flux_error)
#    flux_means = smooth_array(flux_means, 10, 1/spec.flux_error)
    flux_means = smooth_array(flux_means, 5, 1/spec.flux_error)
    flux_means = smooth_array(flux_means, 15, 1/spec.flux_error)
#    flux_means = smooth_array(flux_means, 25, 1/spec.flux_error)

    
    spec_temp = deepcopy(spec)
    print spec_temp.flux
    spec_temp.flux = flux_means
    print spec_temp.flux
    spec_temp.wave = spec.wave
    spec_temp.flux_error = spec.flux_error

    print 'start weird fit'
    fit_powerlaw_individual(spec, settings)
    print 'spec_temp'
    print spec_temp.alpha, spec_temp.delta
    #base_flux = spec_temp.flux / spec_temp.powerlaw
    #base_wave = spec_temp.wave
    base_flux = flux_means
    base_wave = spec.wave

#    print compspec_flux[800:900]

    index = np.where(base_wave)[0]

    # we use 2 emission free regions
    emfree_regions_num = int(np.size(emfree)/2)
    # create emtpy lists for the log wavelength, log flux and error data
    npix = spec.npix
    wave = spec.wave
    flux = spec.flux
    flux_error = spec.flux_error

    wave_log       = np.zeros((emfree_regions_num, npix))
    flux_log       = np.zeros((emfree_regions_num, npix))
    flux_error_log = np.zeros((emfree_regions_num, npix))
    base_flux_log  = np.zeros((emfree_regions_num, npix))
    wave_temp      = np.zeros(npix)
    flux_temp      = np.zeros(npix)
    flux_error_temp= np.zeros(npix)
    base_flux_temp = np.zeros(npix)


    # create a variable, which counts how many regions contain usable data
    # only fit a function if all 4 regions are usable
    emfree_regions_data = 0

    for i in xrange(emfree_regions_num):
        # Find the element in the wavelength array, which is larger / smaller 
        # than the beginning / end of the emission free region
        wave_em_interval_start = find_element_larger_in_arrays(wave, (1.0 + zem)*emfree[i,0], npix)
        # the -1 at the end takes into account, that the function always returns the bigger
        # value.
        # NOTE: Although the C program states to only use the element in the wavelength array that 
        # is the last element smaller than (1+spec.z)*emfree[i,1], it uses the next one. Thus, for
        # now, we neglect the - 1
        wave_em_interval_end = find_element_larger_in_arrays(wave, (1.0 + zem)*emfree[i,1], npix)# - 1
        # define number of elements

        # include wave_em variables!!!
        wave_start = int(wave_em_interval_start)
        wave_end   = int(wave_em_interval_end)
        index = np.where((flux > 0)             &
                         (flux_error > 0)       &
                         (np.isposinf(flux_error) == False) &
                         (spec.mask_comb == 0))[0]
        index = np.extract((wave_start <= index) & (index <= wave_end), index)

        # TODO: call function, which checks flux_log and flux_error_log for very big deviations.
        # if big, then throw out values. return indices of arrays to use
        if np.size(index) == 0:
            continue
        keep_indices = drop_data_from_intervals(wave, flux, flux_error, index, deviation_factor)

        int_length = np.size(keep_indices)
        interval_index = np.arange(int_length)
        # nmed necessary? check if bigger zero, yes, but otherwise no?

        # introduced, because of weird bug, when indices_for_cspec == [] in np.put below
        if np.size(keep_indices) == 0:
            return 0

        # TODO: circumvent usage of temporary arrays?
        wave_temp       = np.log10(wave)
        np.put(wave_log[i], interval_index, wave_temp[keep_indices])

        flux_temp       = np.log10(flux)
        np.put(flux_log[i], interval_index, flux_temp[keep_indices])

        base_flux_temp   = np.log10(base_flux)
        np.put(base_flux_log[i], interval_index, base_flux_temp[keep_indices])

        flux_error_temp = flux_error / (flux * np.log(10))#np.log10(flux_error)
        np.put(flux_error_log[i], interval_index, flux_error_temp[keep_indices])
        
        if np.size(index) > 0:
            emfree_regions_data += 1

    wave_log       = np.ravel(wave_log)
    index          = np.where(wave_log > 0)[0]
    wave_log       = wave_log[index]

    flux_log       = np.ravel(flux_log)[index]
    flux_error_log = np.ravel(flux_error_log)[index]
    base_flux_log = np.ravel(base_flux_log)[index]
    # DEBUGGING
    # if np.isnan(flux_log).any():
    #     #print 'problem!'
    #     #print flux_log
    #     import sys
    #     sys.exit()
    # if np.isnan(wave_log).any():
    #     #print 'problem2!'
    #     #print wave_log
    #     import sys
    #     sys.exit()
    # if np.isnan(flux_error_log).any():
    #     #print 'problem3!'
    #     #print flux_error_log
    #     import sys
    #     sys.exit()
    
    # Fit a linear function to the log data:
    # also fit if we only have 4 of 5 regions for high redshift objects
    if emfree_regions_data >= 1 and np.size(np.isfinite(wave_log)) > 2 and np.size(np.isfinite(flux_log)) > 2:# >= 4:
        # fit linear function func to out log arrays

        # fit a function of powerlaw + scaling*composite spectrum
        def func(x, a, b, c):
            val = a*x + b + c*base_flux_log
            print val
            return val

        # coeff: fitting parameters
        # pcov: covariance matrix, used to retrieve error for alpha.
        # currently we don't have a chi^2 from this function :/
        coeff, pcov = optimize.curve_fit(func, wave_log, flux_log, p0=(spec.beta, spec.delta, 0.5), sigma=flux_error_log, absolute_sigma=True)

        spec.beta = coeff[0]
        spec.alpha = -spec.beta - 2
        spec.delta = coeff[1]
        spec.scaling = coeff[2]
#        spec.chisq = redchisq
        try:
            # C program seems to take variance as error. Normally would take sqrt of pcov.
            spec.alpha_error = float(pcov[0,0])
        except TypeError:
            print "Fitting problem"

#        base_flux = base_flux[citerator:np.size(spec.wave)+citerator]

        npix_cspec = np.size(base_flux)
        end_val    = np.log10(spec.wave[0])+0.0001*npix_cspec

        # calculate the continuum fit:
        scaling = spec.scaling
        spec.powerlaw = 10.0**(coeff[1] + coeff[0]*np.linspace(np.log10(spec.wave[0]), end_val, npix_cspec) + scaling*np.log10(base_flux))

    # if we don't have 2 usable regions, set everything to 0
    else:
        # Currently leave those values at -999. Seems to be same as in C code then.
        print 'afasddf'
        spec.powerlaw = np.zeros(spec.npix)
        # Think about if 0 is a good value (currently checked in build_base)

    if return_data == 0:
        del(wave_log)
        del(flux_log)
        del(flux_error_log)
        #        del(base_flux_log)
    if return_data == 1:
        return wave_log, flux_log, flux_error_log#, base_flux_log

    spec.normalized = spec.flux / spec.powerlaw
    del(flux_temp)
    del(wave_temp)
    del(flux_error_temp)


def fit_powerlaw_plus_compspec(spec, settings, compspec_flux = 0, return_data = 0, deviation_factor = 3.0, zem_index = -999):
# This function fits the powerlaw to the spectrum, by taking the emission free regions
# emfree and fits a linear function to log-log flux- wavelength data
# In contrast to fit_powerlaw() it doesn't take the median of the intervals, but actually
# works on the individual data points in the intervals

# TODO: properly comment function, important since some things not so obvious
    # define emission free regions
#    emfree = np.array([[1280.0,1292.0],[1312.0,1328.0],[1345.0,1365.0],[1440.0,1475.0], [1610,1790]])

    zem     = spec.z
    z_max   = settings.z_max
    z_min   = settings.z_min
    z_delta = settings.z_delta
#    emfree = np.array([[1445.0, 1455.0], [1973.0, 1983.0]])

    fit_powerlaw_individual(spec, settings)
    
    # currently uselessly two regions, because code doesn't work with one
    # however, don't want to restructure, because eventually probably not the smartest
    # thing to take the whole region!
    emfree = np.array([[1205.0, 1400.0], [1400.0, 1800.0]])

    #    hdu = fitsio.FITS('CalcCompFits/DR10_noBALetc_newVarFits.fits')
    hdu = fitsio.FITS('./CalcCompFits/dr9_compspectrum_2.5-3.fits')

#    print compspec_flux[800:900]

    hdu_correct = fitsio.FITS('dr7list_fits_ind.fits')
    ImageData = hdu[0].read()
    compspec_flux = ImageData[3, 0:5763].copy()
#    compspec_flux = compspec_flux - 1.0
    compspec_wave = ImageData[0, 0:5763].copy()

    ImageData_correct = hdu_correct[0].read()
    hdu0_correct_header = hdu_correct[0].read_header()
    npix_correct = hdu0_correct_header['NAXIS1']
    average_trans = np.empty(5763)
    average_trans.fill(1.0)
    coeff0 = hdu0_correct_header['COEFF0']
    coeff1 = hdu0_correct_header['COEFF1']
    wave_correct = 10**(coeff0 + coeff1*np.arange(npix_correct))
    print compspec_wave[0], compspec_wave[-1]
    print wave_correct[0], wave_correct[-1]
    wave_em_interval_start = find_element_larger_in_arrays(compspec_wave, wave_correct[0], 5763)
    wave_em_interval_end = find_element_larger_in_arrays(compspec_wave, wave_correct[-1], 5763)# - 1
    wave_start = int(wave_em_interval_start)
    wave_end   = int(wave_em_interval_end)
    index = np.where(compspec_wave)[0]
    index = np.extract((wave_start <= index) & (index <= wave_end), index)
    temp = 1 - ImageData_correct[0, 0:npix_correct]
    print np.shape(index), np.shape(temp)
    np.put(average_trans, index, temp)
    print 'wave'
#    print index, compspec_wave[index], wave_correct[index]
#    compspec_flux = compspec_flux * average_trans


    compspec_wave = compspec_wave*(1.0 + zem)

    # we use 2 emission free regions
    emfree_regions_num = int(np.size(emfree)/2)
    # create emtpy lists for the log wavelength, log flux and error data
    npix = spec.npix
    wave = spec.wave
    flux = spec.flux
    flux_error = spec.flux_error

    wave_log       = np.zeros((emfree_regions_num, npix))
    flux_log       = np.zeros((emfree_regions_num, npix))
    flux_error_log = np.zeros((emfree_regions_num, npix))
    compspec_flux_log = np.zeros((emfree_regions_num, npix))
    wave_temp      = np.zeros(npix)
    flux_temp      = np.zeros(npix)
    flux_error_temp= np.zeros(npix)
    compspec_flux_temp = np.zeros(npix)


    # create a variable, which counts how many regions contain usable data
    # only fit a function if all 4 regions are usable
    emfree_regions_data = 0

    for i in xrange(emfree_regions_num):
        # Find the element in the wavelength array, which is larger / smaller 
        # than the beginning / end of the emission free region
        wave_em_interval_start = find_element_larger_in_arrays(wave, (1.0 + zem)*emfree[i,0], npix)
        # the -1 at the end takes into account, that the function always returns the bigger
        # value.
        # NOTE: Although the C program states to only use the element in the wavelength array that 
        # is the last element smaller than (1+spec.z)*emfree[i,1], it uses the next one. Thus, for
        # now, we neglect the - 1
        wave_em_interval_end = find_element_larger_in_arrays(wave, (1.0 + zem)*emfree[i,1], npix)# - 1
        # define number of elements

        # include wave_em variables!!!
        wave_start = int(wave_em_interval_start)
        wave_end   = int(wave_em_interval_end)
        index = np.where((flux > 0)             &
                         (flux_error > 0)       &
                         (np.isposinf(flux_error) == False) &
                         (spec.mask_comb == 0))[0]
        index = np.extract((wave_start <= index) & (index <= wave_end), index)

        # TODO: call function, which checks flux_log and flux_error_log for very big deviations.
        # if big, then throw out values. return indices of arrays to use
        if np.size(index) == 0:
            continue
        keep_indices = drop_data_from_intervals(wave, flux, flux_error, index, deviation_factor)

        int_length = np.size(keep_indices)
        interval_index = np.arange(int_length)
        # nmed necessary? check if bigger zero, yes, but otherwise no?

        # introduced, because of weird bug, when indices_for_cspec == [] in np.put below
        if np.size(keep_indices) == 0:
            return 0

        iterator = int(find_element_larger_in_arrays(spec.wave, compspec_wave[0], spec.npix))
        if spec.wave[iterator] > compspec_wave[0] and iterator > 0:
            iterator -= 1
        citerator = int(find_element_larger_in_arrays(compspec_wave, spec.wave[0], spec.npix))
        indices_for_cspec = np.subtract(keep_indices, iterator)
        indices_for_cspec = np.add(indices_for_cspec, citerator)

        # TODO: circumvent usage of temporary arrays?
        wave_temp       = np.log10(wave)
        np.put(wave_log[i], interval_index, wave_temp[keep_indices])

        flux_temp       = np.log10(flux)
        np.put(flux_log[i], interval_index, flux_temp[keep_indices])

        compspec_flux_temp   = np.log10(compspec_flux)
        np.put(compspec_flux_log[i], interval_index, compspec_flux_temp[indices_for_cspec])

        flux_error_temp = flux_error / (flux * np.log(10))#np.log10(flux_error)
        np.put(flux_error_log[i], interval_index, flux_error_temp[keep_indices])
        
        if np.size(index) > 0:
            emfree_regions_data += 1

    wave_log       = np.ravel(wave_log)
    index          = np.where(wave_log > 0)[0]
    wave_log       = wave_log[index]

    flux_log       = np.ravel(flux_log)[index]
    flux_error_log = np.ravel(flux_error_log)[index]
    compspec_flux_log = np.ravel(compspec_flux_log)[index]
    # DEBUGGING
    # if np.isnan(flux_log).any():
    #     #print 'problem!'
    #     #print flux_log
    #     import sys
    #     sys.exit()
    # if np.isnan(wave_log).any():
    #     #print 'problem2!'
    #     #print wave_log
    #     import sys
    #     sys.exit()
    # if np.isnan(flux_error_log).any():
    #     #print 'problem3!'
    #     #print flux_error_log
    #     import sys
    #     sys.exit()
    
    # Fit a linear function to the log data:
    # also fit if we only have 4 of 5 regions for high redshift objects
    if emfree_regions_data >= 1 and np.size(np.isfinite(wave_log)) > 2 and np.size(np.isfinite(flux_log)) > 2:# >= 4:
        # fit linear function func to out log arrays

        # fit a function of powerlaw + scaling*composite spectrum
        def func(x, a, b, c):
            val = a*x + b + c*compspec_flux_log
            print val
            return val

        # coeff: fitting parameters
        # pcov: covariance matrix, used to retrieve error for alpha.
        # currently we don't have a chi^2 from this function :/
        coeff, pcov = optimize.curve_fit(func, wave_log, flux_log, p0=(spec.beta, spec.delta, 0.5), sigma=flux_error_log, absolute_sigma=True)

        spec.beta = coeff[0]
        spec.alpha = -spec.beta - 2
        spec.delta = coeff[1]
        spec.scaling = coeff[2]
#        spec.chisq = redchisq
        try:
            # C program seems to take variance as error. Normally would take sqrt of pcov.
            spec.alpha_error = float(pcov[0,0])
        except TypeError:
            print "Fitting problem"

        # use the fitted coefficients to calculate the powerlaw continuum
        iterator = int(find_element_larger_in_arrays(spec.wave, compspec_wave[0], spec.npix))
        if spec.wave[iterator] > compspec_wave[0] and iterator > 0:
            iterator -= 1
        citerator = int(find_element_larger_in_arrays(compspec_wave, spec.wave[0], spec.npix))
        index = np.where(spec.wave)[0]
        indices_for_cspec = np.subtract(index, iterator)
        indices_for_cspec = np.add(indices_for_cspec, citerator)
#        compspec_flux = compspec_flux[citerator:np.size(spec.wave)+citerator]

        compspec_flux = compspec_flux[indices_for_cspec]
        npix_cspec = np.size(compspec_flux)
        end_val    = np.log10(spec.wave[0])+0.0001*npix_cspec

        # calculate the continuum fit:
        scaling = spec.scaling
        spec.powerlaw = 10.0**(coeff[1] + coeff[0]*np.linspace(np.log10(spec.wave[0]), end_val, npix_cspec) + scaling*np.log10(compspec_flux))

    # if we don't have 2 usable regions, set everything to 0
    else:
        # Currently leave those values at -999. Seems to be same as in C code then.
        print 'afasddf'
        spec.powerlaw = np.zeros(spec.npix)
        # Think about if 0 is a good value (currently checked in build_compspec)

    if return_data == 0:
        del(wave_log)
        del(flux_log)
        del(flux_error_log)
        #        del(compspec_flux_log)
    if return_data == 1:
        return wave_log, flux_log, flux_error_log#, compspec_flux_log

    spec.normalized = spec.flux / spec.powerlaw
    del(flux_temp)
    del(wave_temp)
    del(flux_error_temp)


################################################################################
############################## Powerlaw ########################################
################################################################################
    


#TODO: check influence of np.zeros in code!
def fit_powerlaw_individual(spec, settings, return_data = 0, deviation_factor = 3.0, zem_index = -999, emfree = []):
# This function fits the powerlaw to the spectrum, by taking the emission free regions
# emfree and fits a linear function to log-log flux- wavelength data
# In contrast to fit_powerlaw() it doesn't take the median of the intervals, but actually
# works on the individual data points in the intervals

# TODO: properly comment function, important since some things not so obvious
    # define emission free regions
    #emfree = np.array([[1280.0,1292.0],[1312.0,1328.0],[1345.0,1365.0],[1440.0,1475.0]])#, [1610,1790]])


    zem     = spec.z
    z_max   = settings.z_max
    z_min   = settings.z_min
    z_delta = settings.z_delta

    if emfree == []:
        emfree_intervals = [[1280.0,1292.0],[1312.0,1328.0],[1345.0,1365.0],[1440.0,1475.0], 
                            [1685, 1715], [1730, 1742], [1805, 1837], [2020, 2055], [2190, 2210]]
        # 1. z_max = 3.6, maybe 4.0
        # 2. z_max = 3.6
        # 3. z_max = 3.6
        # 4. z_max = 3.2
        # 5. z_max = 2.8
        # zem_index: 0 == 2.2 - 2.4; 
        #            1 == 2.4 - 2.8; 
        #            2 == 2.8 - 3.2;
        #            3 == 3.2 - 3.6; 
        #            4 == 3.6 - 4.0;
        #            5 == 4.0 - 4.4
        #            ...
        # if zem_index == -999
        if zem_index == -999:
            emfree = np.array(emfree_intervals)
        else:
            emfree_int_tuples= [(emfree_intervals[0], z_max), (emfree_intervals[1], z_max), 
                                (emfree_intervals[2], z_max), (emfree_intervals[3], z_max),
                                (emfree_intervals[4], 4.0),   (emfree_intervals[5], 3.6  ),
                                (emfree_intervals[6], 3.6),   (emfree_intervals[7], 3.2  ),
                                (emfree_intervals[8], 2.8)]
            # construct emfree array from emfree_int_tuples
            emfree = []
            zem_obj = z_min + zem_index*z_delta + z_delta/2
            for i in xrange(len(emfree_int_tuples)):
                if zem_obj < emfree_int_tuples[i][1]:
                    emfree.append(emfree_intervals[i])
    #        print emfree
    #        if zem_index > 4:
    #            emfree = np.asarray(emfree[4])
            emfree = np.asarray(emfree)

    # we use 4 emission free regions
    emfree_regions_num = int(np.size(emfree)/2)
    # create emtpy lists for the log wavelength, log flux and error data
    npix = spec.npix
    wave = spec.wave
    flux = spec.flux
    flux_error = spec.flux_error

    wave_log       = np.zeros((emfree_regions_num, npix))
    flux_log       = np.zeros((emfree_regions_num, npix))
    flux_error_log = np.zeros((emfree_regions_num, npix))
    wave_temp      = np.zeros(npix)
    flux_temp      = np.zeros(npix)
    flux_error_temp= np.zeros(npix)

    # create a variable, which counts how many regions contain usable data
    # only fit a function if all 4 regions are usable
    emfree_regions_data = 0

    for i in xrange(emfree_regions_num):
        # Find the element in the wavelength array, which is larger / smaller 
        # than the beginning / end of the emission free region
        wave_em_interval_start = find_element_larger_in_arrays(wave, (1.0 + zem)*emfree[i,0], npix)
        # the -1 at the end takes into account, that the function always returns the bigger
        # value.
        # NOTE: Although the C program states to only use the element in the wavelength array that 
        # is the last element smaller than (1+spec.z)*emfree[i,1], it uses the next one. Thus, for
        # now, we neglect the - 1
        wave_em_interval_end = find_element_larger_in_arrays(wave, (1.0 + zem)*emfree[i,1], npix)# - 1
        # define number of elements

        # include wave_em variables!!!
        wave_start = int(wave_em_interval_start)
        wave_end   = int(wave_em_interval_end)
        index = np.where(np.logical_and(flux > 0, flux_error > 0, np.isposinf(flux_error) == False))[0]
        index = np.extract((wave_start <= index) & (index <= wave_end), index)

        # TODO: call function, which checks flux_log and flux_error_log for very big deviations.
        # if big, then throw out values. return indices of arrays to use
        if np.size(index) == 0:
            continue
        keep_indices = drop_data_from_intervals(wave, flux, flux_error, index, deviation_factor)

        int_length = np.size(keep_indices)
        interval_index = np.arange(int_length)
        # nmed necessary? check if bigger zero, yes, but otherwise no?

        # TODO: circumvent usage of temporary arrays?
        wave_temp       = np.log10(wave)
        np.put(wave_log[i], interval_index, wave_temp[keep_indices])

        flux_temp       = np.log10(flux)
        np.put(flux_log[i], interval_index, flux_temp[keep_indices])

        flux_error_temp = flux_error / (flux * np.log(10))#np.log10(flux_error)
        np.put(flux_error_log[i], interval_index, flux_error_temp[keep_indices])
        
        if np.size(index) > 0:
            emfree_regions_data += 1

    wave_log       = np.ravel(wave_log)
    index          = np.where(wave_log > 0)[0]
    wave_log       = wave_log[index]

    flux_log       = np.ravel(flux_log)[index]
    flux_error_log = np.ravel(flux_error_log)[index]
    
    # DEBUGGING
    # if np.isnan(flux_log).any():
    #     print 'problem!'
    #     print flux_log
    #     import sys
    #     sys.exit()
    # if np.isnan(wave_log).any():
    #     print 'problem2!'
    #     print wave_log
    #     import sys
    #     sys.exit()
    # if np.isnan(flux_error_log).any():
    #     print 'problem3!'
    #     print flux_error_log
    #     import sys
    #     sys.exit()
    
    # Fit a linear function to the log data:
    # also fit if we only have 4 of 5 regions for high redshift objects
    #def func(x, a, b):
    #    return a*x + b
    if emfree_regions_data >= 4 and np.size(np.isfinite(wave_log)) > 2 and np.size(np.isfinite(flux_log)) > 2:# >= 4:
        # fit linear function func to out log arrays
        # coeff: fitting parameters
        # pcov: covariance matrix, used to retrieve error for alpha.

        #coeff, pcov = optimize.curve_fit(func, wave_log, flux_log, p0=(-2,5), sigma=flux_error_log, absolute_sigma=True)

        print 'testtest'
        coeff, pcov, redchisq = linfit(wave_log, flux_log, sigmay=flux_error_log, cov=True, chisq=True, relsigma=False, residuals=False)
        # assign coefficients to our spectrum
        spec.beta = coeff[0]
        print coeff[0], coeff[1]
        spec.alpha = -spec.beta - 2
        spec.delta = coeff[1]
        spec.chisq = redchisq
        try:
            # C program seems to take variance as error. Normally would take sqrt of pcov.
            spec.alpha_error = float(pcov[0,0])
        except TypeError:
            print "Fitting problem"

        # use the fitted coefficients to calculate the powerlaw continuum
        spec.powerlaw = 10.0**(coeff[1] + coeff[0]*np.log10(spec.wave))

    # if we don't have 4 usable regions, set everything to 0
    else:
        # Currently leave those values at -999. Seems to be same as in C code then.
        print 'ok?'
        spec.powerlaw = np.zeros(spec.npix)
        # Think about if 0 is a good value (currently checked in build_compspec)

    if return_data == 0:
        del(wave_log)
        del(flux_log)
        del(flux_error_log)
    if return_data == 1:
        return wave_log, flux_log, flux_error_log


    spec.normalized = spec.flux / spec.powerlaw
    del(flux_temp)
    del(wave_temp)
    del(flux_error_temp)

def fit_powerlaw(spec, return_data = 0):
# This function fits the powerlaw to the spectrum, by taking the emission free regions
# emfree and fits a linear function to log-log flux- wavelength data
    emfree = []
    # define emission free regions
    emfree = np.array([[1280.0,1292.0],[1312.0,1328.0],[1345.0,1365.0],[1440.0,1475.0]])#, [1610,1790]])
    # we use 4 emission free regions
    emfree_regions_num = 4#5
    # create emtpy lists for the log wavelength, log flux and error data
    wave_log       = np.zeros(emfree_regions_num)
    flux_log       = np.zeros(emfree_regions_num)
    flux_error_log = np.zeros(emfree_regions_num)
    flux_temp      = np.zeros(np.size(spec.flux))
#    flux_error_temp= np.zeros(np.size(spec.flux_error))

    # set standard values for alpha, alpha_error, beta and delta
    spec.alpha       = -999
    spec.alpha_error = -999
    spec.beta        = -999
    spec.delta       = -999
    spec.chisq       = -999
    # create a variable, which counts how many regions contain usable data
    # only fit a function if all 4 regions are usable
    emfree_regions_data = 0

    # Initialise spectrum's emfree matrix:
    spec.emfree = np.zeros((3,emfree_regions_num))
    for i in xrange(emfree_regions_num):
        # calculate the array containing the emission free regions
        # TODO: take spec.emfree calculations out of for loop. Probably hardly speed difference
        spec.emfree[0,i] = (1.0 + spec.z)*0.5*(emfree[i,0] + emfree[i,1])
        spec.emfree[1,i] = 0.0
        spec.emfree[2,i] = -1.0

        # Find the element in the wavelength array, which is larger / smaller 
        # than the beginning / end of the emission free region
        wave_em_interval_start = find_element_larger_in_arrays(spec.wave, (1.0 + spec.z)*emfree[i,0], spec.npix)
        # the -1 at the end takes into account, that the function always returns the bigger
        # value.
        # TODO: Check why siqr not the same as in C code and median neither. (slight abbreviations)
        # NOTE: Although the C program states to only use the element in the wavelength array that 
        # is the last element smaller than (1+spec.z)*emfree[i,1], it uses the next one. Thus, for
        # now, we neglect the - 1
        wave_em_interval_end = find_element_larger_in_arrays(spec.wave, (1.0 + spec.z)*emfree[i,1], spec.npix)# - 1
        # define number of elements
        # nmed is the number of usable pixels in the flux array. It is given by the number of pixels
        # in the interval of wave_em_interval_start to end with a flux error value bigger than 0
        nmed = 0
        for k in xrange(wave_em_interval_start, wave_em_interval_end):
            # if there is a non-vanishing error for the pixel, it means that
            # the pixel is good to use. Save all those flux elements for calc
            # of median in temporary flux array and count number of pixels, nmed
            if spec.flux_error[k] > 0:
                flux_temp[nmed]       = spec.flux[k]
                nmed += 1

        # Calculation of semi-interquartile using built in numpy functions. I can't manage to get the
        # same results as the C program by using them, so I implemented my own function
        # The reason is probably how the elements are chosen. 
        # in numpu 1.8 we can't choose different interpolations. Possible from 1.9
        # percentile84 = np.percentile(spec.flux[wave_em_interval_start:wave_em_interval_end], 84)
        # percentile16 = np.percentile(spec.flux[wave_em_interval_start:wave_em_interval_end], 16)
        # siqr = (percentile84 - percentile16)/2.0
        
        if nmed > 0:
            median = np.median(flux_temp[0:nmed])
            siqr = calc_siqr(flux_temp[0:nmed], nmed)
            
            spec.emfree[1,i] = median
            spec.emfree[2,i] = siqr/np.sqrt(nmed)
            # if usable values, append emfree regions to log data arrays
            if spec.emfree[1,i] > 0 and spec.emfree[2,i] > 0:
                wave_log[i] = np.log10(spec.emfree[0,i])
                flux_log[i] = np.log10(spec.emfree[1,i])
                flux_error_log[i] = np.log10(1.0 + spec.emfree[2,i]/spec.emfree[1,i])
                # count region as usable
                emfree_regions_data += 1
    
    # Fit a linear function to the log data:
    # define our (line) fitting function
    # using linfit we don't need our linear fitting function
    # Linfit gives same results as optimize.curve now! 
    # def func(x, a, b):
    #    return a*x + b
#    only continue, if 4 regions contain usable data
    if emfree_regions_data == emfree_regions_num:# >= 4:
        # fit linear function func to out log arrays
        # coeff: fitting parameters
        # pcov: covariance matrix, used to retrieve error for alpha.
#        coeff, pcov = optimize.curve_fit(func, wave_log, flux_log, p0=(-2,5), sigma=flux_error_log)

#        print 'normale func:', wave_log, flux_log

        coeff, pcov, redchisq, residuals = linfit(wave_log, flux_log, sigmay=flux_error_log, cov=True, chisq=True, relsigma=False, residuals=True)
        
#        print wave_log, flux_log
        
        # assign coefficients to our spectrum
        spec.beta = coeff[0]
        spec.alpha = -spec.beta - 2
        spec.delta = coeff[1]
        spec.chisq = redchisq
        try:
            # C program seems to take variance as error. Normally would take sqrt of pcov.
            spec.alpha_error = float(pcov[0,0])
        except TypeError:
            print "Fitting problem"

#        print "alpha, delta: ", spec.alpha, spec.delta, spec.alpha_error
        # use the fitted coefficients to calculate the powerlaw continuum
        spec.powerlaw = 10.0**(coeff[1] + coeff[0]*np.log10(spec.wave))

    # if we don't have 4 usable regions, set everything to 0
    else:
        # Currently leave those values at -999. Seems to be same as in C code then.
        # spec.beta = 0
        # spec.alpha = 0
        # spec.delta = 0
        # spec.alpha_error = 0
        spec.powerlaw = np.zeros(spec.npix)
        # Think about if 0 is a good value (currently checked in build_compspec)

    if return_data == 0:
        del(wave_log)
        del(flux_log)
        del(flux_error_log)
    if return_data == 1:
        return wave_log, flux_log, flux_error_log

    del(spec.emfree)
    del(flux_temp)
    # TODO: values differ slightly from c program!
    # Fixed by implementing own siqr function. fixed again, alpha this time
