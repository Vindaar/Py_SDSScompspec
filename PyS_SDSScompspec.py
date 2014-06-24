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
from guppy import hpy

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
    
    if args_check(args, settings) == 0:
        return 0


    # Basic declarations
    files = list(settings.inputfile)
    dustmap = '/home/basti/SDSS_indie/dust_maps/maps/SFD_dust_4096_%s.fits'
    spectra = np.array([spectrum() for i in xrange(len(files))])
    compspec = comp_spectrum(5763)
    i = 0
    a = create_colorcurves()
    alpha_top = 1.5
    alpha_low = -2
    spectra_count = 0
    # memory analysis:

    # Read filename for the output FITS file:
    if settings.outfile == '':
        settings.outfile = raw_input('Give the name of the output FITS file: ')

    # Start the loop over all files in the 
    for i, file in enumerate(files):
        print "Starting with spectrum #: ", i
        filetype = check_filetype(file)
        if filetype == 1:
            # If the return value is one, we do a check on coordinates and this 
            # object is outside the wanted boundaries. Therefore neglect and continue
            if read_spSpec(file, spectra[i], settings) == 1:
                continue
        if filetype == 2:
            if read_spec(file, spectra[i], settings) == 1:
                continue
        
        # Conditions on the QSOs:
        if spectra[i].z > 2.2 and spectra[i].z < 5.3:
            # Dust corrections. Only done, if --dust flag is set on startup
            if settings.dust == 1:
                get_Ebv(spectra[i], dustmap)
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
        # Every 50th loop, we free all objects, which are not used anymore. 
        # This function is called automatically, but not often enough. Reduces
        # memory usage quite a lot.
        if i % 50 == 0:
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
    build_fits_file(compspec, spectra, settings.outfile, settings)
    print "Spectra used: ", compspec.spectra_count, "/", len(files)

    # If cspec flag is set, we return the compspec object to the program that called this function
    if settings.cspec:
        return compspec

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])


################################# END OF THE PROGRAM #################################
