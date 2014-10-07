# this module contains general utility functions used in the whole project,
# like command line argumetn checks, help screen
# as well as computational utility functions
# like align two arrays etc.

import numpy as np

# re contains the function search, which is used to search through
# argv
import re

#####################################################################################
################################### Array handling ##################################
#####################################################################################


def smooth_array(array, smooth_width, weights = None):
    nsize = int(np.size(array))
    # if weights == None:
    #     weights = np.empty(nsize)
    #     weights.fill(1.0)
    array_temp = np.zeros(nsize)
    array_orig = deepcopy(array)
    if weights != None:
        mean_snr = np.ma.average(array * weights)
        snr = array * weights
        #array = array * (array * weights)
        if mean_snr > 10.0:
            upper = 5.0
            lower = 0.3
        else:
            upper = mean_snr/2.7
            lower = 3.0*(np.tanh(mean_snr) - 1.0)
        
        
    for i in xrange(np.size(array)):
        start = end = 0
        if i - smooth_width < 0:
            start = 0
        else:
            start = i - smooth_width
        if i + smooth_width > nsize:
            end = nsize
        else:
            end   = i + smooth_width
        if weights != None:
            std = np.std(array[start:end])
            mean = np.ma.average(array[start:end])#, weights=weights[start:end])
            index = np.where(np.absolute(array[start:end] - mean) < 2.0*std)[0]
            #index = np.where((array[start:end] - mean < 5.0*std) & 
                             #(array[start:end] - mean > -2.3*std))[0]
                             #(array[start:end] - mean > -ratio*std         ) &
                             #(array[start:end] - mean < -ratio/2.0*std         ))&[0]
            index = np.add(index, start)
            array_temp[i] = np.ma.average(array_orig[index], weights=1/snr[index])
#            array_temp[i] = array_temp[i] / (array_orig[i] * weights[i])
        else:
            std = np.std(array[start:end])
            mean = np.average(array[start:end])
            index = np.where(np.absolute(array[start:end] - mean) < 0.5*std)[0]
            index = np.add(index, start)
            array_temp[i] = np.average(array[index])

    if weights != None:
#        array_temp = array_temp / mean_snr
        array      = array_orig
    return array_temp

def find_element_larger_in_arrays(wave, target_wave, npix):
    # Function which takes in a wavelength array wave and a target wavelength
    # It searches for the i-th element in wave, which is the first element
    # bigger than target_wave. 
    # Note: i is always larger than target_wave! If smaller than target_wave is 
    # wanted, take returnvalue - 1

    # Numpy where MUCH faster than using while!

    if wave[npix-1] < target_wave:
        return npix-1
    else:
        i = np.where(wave > target_wave)[0][0]
        return i

    # old depcrecated way of calculating find_element
    # while True:
    #     if i < npix and wave[i] < target_wave:
    #         i += 1
    #     else:
    #         break


def calc_siqr(flux, nmed):
    # Function which calculates the 68% semi-interquartile range, in the exact same way
    # as done in the C program.
    # range_val is the range of the 68% semi-interquartile
    range_val = int(0.341344746069*nmed)
    # we sort the flux array
    flux = np.sort(flux)
    # and define the rounded half of the number of flux elements
    sidx = int(nmed/2)
    if nmed%2 == 0:
        siqr = 0.25*(-flux[sidx - range_val - 1] - flux[sidx - range_val] + flux[sidx+range_val - 1] + flux[sidx+range_val])
    elif nmed%2:
        siqr = 0.5*(flux[sidx+range_val] - flux[sidx - range_val])
    return siqr

def drop_data_from_intervals(wave, flux, flux_error, index, deviation_factor):
    # working on log data currently
    if np.size(index) > 0:
        try:
            mean       = np.average(flux[index], weights=(1/flux_error[index]))
        except ZeroDivisionError:
            # currently we just drop the whole interval
            return []
    else:
        return []
#    mean_error = np.mean(flux_error[index])
    std        = np.absolute(np.std(flux[index]))
#    std_error  = np.absolute(np.std(flux_error[index]))

    # deviation_factor is the factor that the flux allows is allowed to deviate from
    # the mean in the interval
    indices = np.where(np.absolute(flux - mean) < deviation_factor*std)[0]
#    indices = np.where(np.logical_and(np.absolute(flux - mean) < 3*std, np.absolute(flux_error - mean_error) < 3*std_error))[0]
    indices = np.intersect1d(index, indices)
    return indices

def calc_zem_index(zem, z_min, delta):
    # index for binning data into redshift bins
    if zem != 999:
        zem_index = (zem - z_min)/delta
    else:
        zem_index = -1
    if zem_index > 0:
        zem_index = int(zem_index)
    return zem_index




################################################################################
############################## Input / Argschecks ##############################
################################################################################

def args_check(args, settings):
    # Function which checks for command line arguments
    
    if len(args) > 0:
        if (re.search('.fits', args[0]) or re.search('.fit', args[0])):
            print "Working on single FITS files currently not supported."            
            print "Please provide a file containing a list of FITS files."
            return 0
        else:
            # TODO: include exception handling for wrong files?
            settings.inputfile = open(args[0], 'r')
    else:
        help()
        return 0

    if '--dust' in args:
        settings.dust = 1
    else:
        settings.dust = 0
    if '--coords' in args:
        print args
        i = args.index('--coords')
        print i, args
        try:
            settings.l_min  = float(args[i+1])
            settings.b_min  = float(args[i+2])
            settings.l_max  = float(args[i+3])
            settings.b_max  = float(args[i+4])
            settings.coords = 1
        except IndexError:
            print 'Error: Full set of coordinates have to be provided in the following way:'
            print 'l_min b_min l_max b_max'

    if ('-h' or '--help') in args:
        help()
        return 0
    if '-o' in args:
        i = args.index('o')
        try:
            settings.outfile = args[i+1]
        except IndexError:
            print "Error: if -o, --output is set, need to provide an output file name!"
            return 0
    if '--output' in args:
        i = args.index('--output')
        try:
            settings.outfile = args[i+1]
        except IndexError:
            print "Error: if -o, --output is set, need to provide an output file name!"
            return 0
    if '--cspec' in args:
        settings.cspec = 1
    if '--nprocs' in args:
        i = args.index('--nprocs')
        try:
            settings.nprocs = eval(args[i+1])
        except IndexError:
            print "Error: if --nprocs is set, need to provide number of processes to use!"
            return 0
    if '--deviaton_factor' in args:
        i = args.index('--deviation_factor')
        try:
            settings.deviation_factor = eval(args[i+1])
        except IndexError:
            print "Error: if --deviation_factor is set, need to provide a value!"
            return 0
    if '--nside' in args:
        i = args.index('--nside')
        try:
            settings.map_nside = eval(args[i+1])
            settings.map_npix  = 0
        except IndexError:
            print "Error: if --nside is set, need to provide number of processes to use!"
            return 0
    if '--delta_l' in args:
        i = args.index('--delta_l')
        try:
            settings.delta_l = eval(args[i+1])
            if '--delta_b' not in args:
                settings.delta_b / settings.delta_l / 2.0
        except IndexError:
            print "Error: --delta_l set, but no value given!"
            return 0
    if '--delta_b' in args:
        i = args.index('--delta_b')
        try:
            settings.delta_b = eval(args[i+1])
        except IndexError:
            print "Error: --delta_b set, but no value given!"
            return 0
    if '--flux_corr' in args:
        i = args.index('--flux_corr')
        try:
            settings.flux_corr_list = open(args[i+1], 'r').readlines()
            settings.flux_corr = 1
        except IndexError:
            print "Error: if you set --flux_corr flag, you have to provide"
            print "a list of FITS files, which contain the flux correction"
            print "functions."
            return 0


################################################################################
############################## Help screen #####################################
################################################################################


def help():
    print "####################    SDSScompspec in Python    ####################"
    print ""
    print "Usage:"
    print "./PyS_SDSScompspec <input list file> --options"
    print "possible options:"
    print "    --dust:       activates galactic dust corrections"
    print "    -o, --output: provide name for output FITS file"
    print "    --coords:     input a cut on galactic coordinates in degrees"
    print "                  provide in the following way:"
    print "                  l_min b_min l_max b_max"
    print "    --cspec:      activates return of compspec object at the end" 
    print "                  of the program."  
    print "    --deviation_factor: "
    print "                  factor, which determines num*std flux is allowed to"
    print "                  deviate from mean flux in continuum fit intervals"
    print "    --delta_l:    size of considered rectangles in SDSScoordinates in"
    print "                  in latitude"
    print "    --delta_b:    size of considered rectangles in SDSScoordinates in"
    print "                  in longitude. If --delta_l set, but --delta_b not,"
    print "                  we assume delta_b = delta_l / 2"
    print "    -h, --help:   print this help"


