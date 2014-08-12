#!/usr/bin/env python

# Import used modules
# numpy: used for data storage etc.
import numpy as np
# gc contains gc.collect(), which is the function to collect all 
# unused memory (called automatically, just want to call it slightly
# more often
import gc
# import sys
# from matplotlib import rc
# import math
# # the multiprocessing module
import multiprocessing as mp
# guppy to track memory usage. Currently not used
# from guppy import hpy


from SDSSmodules import *
from SDSSclasses import *

# HDU numbers in this file are given starting from 0, since it is
# conform with astropy, although the FITS standard starts from 1!


#####################################################################################    
############################### RUN THE PROGRAM #####################################
#####################################################################################

# TODO: put argcheck outside of function into SDSSmodule so that I don't have to do
# it for both programs.

def main(args, settings = program_settings()):
    # Create settings object, which store flags, which are set on program start
#    settings = program_settings()
    # this way the settings argument should be optional. If we provide it in PyS_SDSScoordinates
    # it should not be necessary here.

    #TODO: IF we supply settings object. Then maybe check if we give input file
    # IF we do that, simply bypass all args checks, since the assume that everything
    # is already given in settings object.

    # Read in filename and settings from command line
    print 'Checking command line arguments'
    if args_check(args, settings) == 0:
        return 0

    # Basic declarations
    if settings.spectra_list:
        files = settings.spectra_list
    else:
        files = list(settings.inputfile)
    dustmap = '/home/basti/SDSS_indie/dust_maps/maps/SFD_dust_4096_%s.fits'
    print 'creating array of spectrum objects'
    nspec = len(files)
    spectra = np.array([spectrum() for i in xrange(nspec)])
    print 'creating object storing coordinate arrays'
    coordinates = coordinate_arrays(nspec)

    # Create a compspec object with 5763 pixels. Taken from C code
    compspec = comp_spectrum(5763)
    i = 0
    print 'create color curves'
    a = create_colorcurves()
    alpha_top = 1.5
    alpha_low = -2
    spectra_count = 0
    # memory analysis:

    files_used_file = open('files_used', 'w')
    files_not_used  = open('files_not_used.txt', 'w')
    files_used = []
    files_not_used_zem = []
    files_not_used_alpha = []
    alpha_of_files_not_used = []
    files_not_used_compspec = []

    alpha_wrong_count = 0

    # Read filename for the output FITS file:
    if settings.outfile == '':
        settings.outfile = raw_input('Give the name of the output FITS file: ')

    # Do filetype check only once
    filetype = check_filetype(files[0])

    # Start the loop over all files in the 
    for i, file in enumerate(files):
        if i % 5000 == 0:
            # this buffer is used to read a bunch of fits files after another,
            # because that way the reading is much faster. Needs more memory
            # of course.
            if i+5000 < len(files):
                buffer = 5000
            else:
                buffer = len(files) - i
            for j in xrange(buffer):
                if j % 100 == 0:
                    print "Reading spectrum #: ", j
                if filetype == 1:
                    # If the return value is one, we do a check on coordinates and this 
                    # object is outside the wanted boundaries. Therefore neglect and continue
                    if read_spSpec_fitsio(files[i+j], spectra[i+j], settings) == 1:
                        continue
                if filetype == 2:
                    if read_spec_fitsio(files[i+j], spectra[i+j], settings) == 1:
                        continue
            # Fill the l and b arrays with the coordinates of the buffer
            # and get the E(B-V) values from the dustmap
            print 'filling coordinate arrays from buffer...'
            fill_coordinate_arrays_from_buffer(coordinates, spectra, dustmap, i, buffer)
        if i % 100 == 0:
            print "Working on spectrum #: ", i
        spectra[i].filename = file
        # Conditions on the QSOs:
        if spectra[i].z >= 2.2 and spectra[i].z <= 5.3:
            # Dust corrections. Only done, if --dust flag is set on startup
            if settings.dust == 1:
                Gal_extinction_correction(spectra[i])

            # Build the median array? Use: flux, continuum, npix,
            # status
            colors(spectra[i], a)
            #powerlaw function
            # TODO: check in powerlaw if QSO usable?
            fit_powerlaw(spectra[i])
            # alpha cut
            if spectra[i].alpha < alpha_top and spectra[i].alpha > alpha_low:
                # calculate comp spec
                if build_compspec(compspec, spectra[i]) == 0:
                    compspec.spectra_count += 1
                    # flag spectrum as used and contributed
                    spectra[i].flag = 0
                    files_used.append(file)
                else:
                    files_not_used_compspec.append(file)
                    # flag spectrum as not used, discarded in build_compspec
                    spectra[i].flag = 3
            else:
                files_not_used_alpha.append(file)
                alpha_of_files_not_used.append(spectra[i].alpha)
                # flag spectrum as not used, because of cut on spectral index
                spectra[i].flag = 2
        else:
            files_not_used_zem.append(file)
            # flag spectrum as not used, because of redshift
            spectra[i].flag = 1
        # Every 50th loop, we free all objects, which are not used anymore. 
        # This function is called automatically, but not often enough. Reduces
        # memory usage quite a lot.
        if spectra[i].alpha == -999:
            # print ""
            # print ""
            # print file, spectra[i].alpha, spectra[i].z
            alpha_wrong_count += 1
        if i % 500 == 0:
            gc.collect()
        # Free all big arrays, which won't be needed anymore, after this loop. 
        # Unecessary memory usage.
        del(spectra[i].flux)
        del(spectra[i].flux_error)
        del(spectra[i].wave)
        del(spectra[i].powerlaw)
        del(spectra[i].status)
        del(spectra[i].snr)
        # increase loop variable, since we run over files and not an integer
        i += 1
	        
    statistics(compspec, spectra)
    # TODO: change build_fits_file such that it reads outfile from settings object within function
    if settings.cspec == 0:
        build_fits_file(compspec, spectra, settings.outfile, settings)

    for i in xrange(len(files_used)):
        files_used_file.write(files_used[i])
    print "Spectra used: ", compspec.spectra_count, "/", len(files)

    print "There are ", alpha_wrong_count, "objects with alpha -999"

    for item in files_not_used_zem:
        print>>files_not_used, "zem     ", item
    for i, item in enumerate(files_not_used_alpha):
        print>>files_not_used, "alpha     ", alpha_of_files_not_used[i], item
    for item in files_not_used_compspec:
        print>>files_not_used, "compspec     ", item

    # If cspec flag is set, we return the compspec object to the program that called this function
    if settings.cspec:
        return compspec

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])


################################# END OF THE PROGRAM #################################
