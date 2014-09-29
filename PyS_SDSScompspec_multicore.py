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
# Guppy to track memory usage. Currently not used
#from guppy import hpy

from SDSSmodules import *
from SDSSclasses import *

# HDU numbers in this file are given starting from 0, since it is
# conform with astropy, although the FITS standard starts from 1!


# TODO:!!!
# in the output file (ran on dr7small.lis) there are 5 objects with 0 everywhere.
# check why and if always 5 or percentage


#####################################################################################    
############################### RUN THE PROGRAM #####################################
#####################################################################################

def file_loop(files, compspec, start_iter, end_iter, out_q, out_q2, out_qcoord_temp, settings):
    # some basic variable definitions and values, which are used to cut on
    # alpha
    i = 0
    j = 0
    dustmap = '/home/basti/SDSS_indie/dust_maps/maps/SFD_dust_4096_%s.fits'
    # a is the array, which contains the color curves. This way it is called for 
    # each process once. Could be moved outside into main function so it is
    # only called once in total, but then we would have to pass a copy to
    # each process.
    a = create_colorcurves()
    alpha_top = 1.5
    alpha_low = -2
    nspec = len(files)

    # start_iter is the first element a process will work on
    i = start_iter
    # create a numpy array of spectrum objects
    spectra = np.array([spectrum() for i in xrange(end_iter - start_iter)])

    # the loop will run over all files from start_iter to end_iter
    for i in xrange(start_iter, nspec):
        if i < end_iter:
            print "Starting with spectrum #: ", i
            # Check for filetype (spSpec and spec) and call appropriate function
            filetype = check_filetype(files[i])
            if filetype == 1:
                read_spSpec(files[i], spectra[j], settings)
            if filetype == 2:
                read_spec_fitsio(files[i], spectra[j], settings)
            elif filetype == 0:
                print "File ", files[i], " neither spSpec nor spec FITS file"
                print "Continue with next file"
                continue

            # Conditions on the QSOs:
            if spectra[j].z > 2.2 and spectra[j].z < 5.3:
                # Start the dust correction. Only done if --dust flag is set
                if settings.dust == 1:
                    get_Ebv(spectra[j], dustmap)
                    Gal_extinction_correction(spectra[j])
                # Build the median array? Use: flux, continuum, npix, status?
                colors(spectra[j], a)
                # powerlaw function
                fit_powerlaw(spectra[j])
                # alpha cut
                if spectra[j].alpha < alpha_top and spectra[j].alpha > alpha_low:
                    # add current spectrum to compspec object
                    if build_compspec(compspec, spectra[j]) == 0:
                        # Count the used QSOs, if QSO was used.
                        compspec.spectra_count += 1
                        # flag spectrum as used and contributed
                        spectra[j].flag = 0
                    else:
                        # flag spectrum as not used, discarded in build_compspec
                        spectra[j].flag = 3
                else:
                    # flag spectrum as not used, because of cut on spectral index
                    spectra[j].flag = 2
            else:
                # flag spectrum as not used, because of redshift
                spectra[j].flag = 1
            # Every 50th loop, we free all objects, which are not used anymore. 
            # This function is called automatically, but not often enough. Reduces
            # memory usage quite a lot.
            if j % 50 == 0:
                gc.collect()
            # Free all big arrays, which won't be needed anymore, after this loop. 
            # Unecessary memory usage.
            del(spectra[j].flux)
            del(spectra[j].flux_error)
            del(spectra[j].wave)
            del(spectra[j].powerlaw)
            del(spectra[j].status)
            del(spectra[j].snr)
#            spectra[j].coordinates = 0
            i += 1
            j += 1
        else: 
            break
    # Add the dictionary of spectra and the compspec to two queues, which will 
    # eventually be combined from all processes

    # Workaround for now? Return a list at the same time containing the coordinates and
    # reassign those coordinates in main function?
    # two dimensional array
    # 0: ra.deg, 1: dec.deg
    # TODO: temporary hotfix, because SkyCoord fails to return self on call of
    # out_q.get(). 
    coordinate_temp_array = np.zeros((len(spectra), 2))
    
    for i in xrange(len(spectra)):
        # coordinates are in FK5
        coordinate_temp_array[i,0] = spectra[i].coordinates.ra.deg
        coordinate_temp_array[i,1] = spectra[i].coordinates.dec.deg
        spectra[i].coordinates = 0

    out_q.put(spectra)
    out_q2.put(compspec)
    out_qcoord_temp.put(coordinate_temp_array)

def main(args, settings = program_settings()):
    # Create settings object, which saves settings set on start
    settings.program_name = 'PyS_SDSScompspec_multicore'

# Read in filename and settings from command line
    if args_check(args, settings) == 0:
        return 0

    print "The program will be run with", settings.nprocs, "subprocesses."
    # Basic declarations
    # Create a list of all files contained in input file; easier to work with

    if settings.spectra_list:
        files = settings.spectra_list
    else:
        files = list(settings.inputfile)
    nspec = len(files)
    # Create a compspec object with 5763 pixels
    compspec = comp_spectrum(5763)
    # Read filename for the output FITS file:
    if settings.outfile == '':
        settings.outfile = raw_input('Give the name of the output FITS file: ')

    # Two Queues are created, which are used for parallel processing. 
    # They will contain:
    # out_q:  dictionary of spectra, returned by each process
    # out_q2: a queue of compspectra (one from each process), which will
    #         be combined by the combine_compspecs() function at the end
    # out_qcoord_temp: a queue containing an array for each process, which
    # contains the coordinates of the corresponding spectra. index i of spectrum
    # list corresponds to index i in coord_array.
    # TODO: this is a temporary hotfix, because SkyCoord object fails to return
    # self on out_q.get() before joining processes.
    out_q           = mp.Queue()
    out_q2          = mp.Queue()
    out_qcoord_temp = mp.Queue()
    # an empty list containing all processes
    procs = []
    # Calculate the number of files each process has to work on
    chunksize = int(math.ceil(nspec / float(settings.nprocs)))
    print nspec, chunksize
    
    # loop over the number of processes and create a new process

    print 'Start parallel processes'
    for i in xrange(settings.nprocs):
        # create a new process
        p = mp.Process(target=file_loop, 
                       args=(files, compspec, chunksize*i, chunksize*(i+1), out_q, out_q2, out_qcoord_temp, settings))
        # add it to the list of processes
        procs.append(p)
        # start the process
        p.start()

    # Create an empty dictionary for the spectra and an empty list for 
    # the compspectra
    # TODO: use numpy arrays instead
    result_spectra = []
    result_compspec = []
    # run over the number of processes in order to generate a final
    # dictionary containing all spectra and a list of all compspecs
    print 'Start combining results of parallel processes'
    for i in xrange(settings.nprocs):
        try:
            result_spectra.append(out_q.get())
        except Empty:
            print i, 'Empty exception raised on element'
        # get next object from out_qcoord_temp queue and save those
        # coordinates in spectra.coordinate object (SkyCoord)
        coordinate_temp_array = out_qcoord_temp.get()
        # check length of array, run over all elements and assign coordinate
        # object to spectrum
        for j in xrange(len(coordinate_temp_array)):
            ra  = coordinate_temp_array[j,0]*u.degree
            dec = coordinate_temp_array[j,1]*u.degree
            # i: access i-th numpy array
            # j: access j-th element in the numpy-array -> a spectrum
            result_spectra[i][j].coordinates = SkyCoord(ra = ra, dec = dec, frame='fk5')
        result_compspec.append(out_q2.get())
    # wait for all processes to finish
    print 'Finished appending items from queues to lists, starting to join...'
    for p in procs:
        p.join()

    print 'Finished joining processes'
    # result_spectra now contains a list of nprocs numpy arrays, which contains the spectra.
    # We want a numpy array of spectra instread:
    # TODO: Think of nicer way to do this
    result_spectra = np.resize(result_spectra, (nspec))

    # calculate the combined compspec from the processes individual ones
    compspec_sum = combine_compspecs(result_compspec)
    # calculate some statistics, which will be printed to the output FITS file
    statistics(compspec_sum, result_spectra)
    # Create the ouput FITS file
    if settings.cspec == 0:
        build_fits_file(compspec_sum, result_spectra, settings.outfile, settings)
    print "Spectra used: ", compspec_sum.spectra_count, "/", nspec

    # If cspec flag is set, we return the compspec object to the program that called this function
    if settings.cspec:
        return compspec_sum

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])


################################# END OF THE PROGRAM #################################
