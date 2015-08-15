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


from SDSSmodules.SDSSfiles import *
from SDSSmodules.SDSSfitting import *
from SDSSmodules.SDSScorrections import *
from SDSSmodules.SDSSoutput import *
from SDSSmodules.SDSScspec import *
from SDSSmodules.SDSSutilities import *
from SDSSmodules.SDSSclasses import *

import cPickle
import theano
import theano.tensor as T
from handle_data import load_SDSS_predict
from mlp import MLP
from convolutional_mlp import LeNetConvNetwork
from congrid import rebin_1d


# HDU numbers in this file are given starting from 0, since it is
# conform with astropy, although the FITS standard starts from 1!


#####################################################################################    
############################### RUN THE PROGRAM #####################################
#####################################################################################

def main(args, settings = program_settings()):
    # Create settings object, which store flags, which are set on program start
#    settings = program_settings()
    # this way the settings argument should be optional. If we provide it in PyS_SDSScoordinates
    # it should not be necessary here.

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

    # Read filename for the output FITS file:
    if settings.outfile == '':
        settings.outfile = raw_input('Give the name of the output FITS file: ')

    nspec = len(files)
#    spectra = np.array([spectrum() for i in xrange(nspec)])

    # Create a compspec object with 5763 pixels. Taken from C code
    compspec = comp_spectrum(5763)
    i = 0
    #print 'create color curves'
    #a = create_colorcurves()
    alpha_top = 1.5
    alpha_low = -2
    spectra_count = 0
    z_min     = settings.z_min
    z_max     = settings.z_max
    z_delta   = settings.z_delta

    # memory analysis:

    files_used_file = open('files_used', 'w')
    files_not_used  = open('files_not_used.txt', 'w')
    files_used = []
    files_not_used_zem = []
    files_not_used_alpha = []
    alpha_of_files_not_used = []
    files_not_used_compspec = []

    # to compare alpha of real and artificial data:
    # simple list that's going to be pickled at the end.
    alpha_of_files_used = []


    alpha_wrong_count = 0

    # resid correction:
    resid_corr = read_resid_corr('residcorr_v5_4_45.dat')

    # Do filetype check only once
    filetype = check_filetype(files[0])

    # create necessary stuff for ANN

    x = T.matrix('x')
    y = T.matrix('y')
    index = T.lscalar()

    rng = np.random.RandomState(1234)

    # read ANN parameters from file:
    if settings.predict_from_mlp == True:
        saved_ann = open(settings.saved_ann)
        layer_params = cPickle.load(saved_ann)
        ann_layout   = cPickle.load(saved_ann)
        saved_ann.close()
            
    # construct the MLP class
    if '--mlp' in args:
        classifier = MLP(
            rng=rng,
            input=x,
            layer_params=layer_params,
            ann_layout=ann_layout,
            activation=T.tanh
        )
    
    # construct the LeNet class
    if '--lenet' in args:
        nkerns = [20, 50]
        batch_size = 1
        #layer0_input = x.reshape((batch_size, 1, 24, 24))
        layer0_input = x.reshape((batch_size, 1, 44, 44))
        
        classifier = LeNetConvNetwork(
            rng=rng,
            input=layer0_input,
            filter0_shape=(nkerns[0], 1, 5, 5),
            #image0_shape=(batch_size, 1, 24, 24),
            image0_shape=(batch_size, 1, 44, 44),
            poolsize=(2,2),
            nkerns=nkerns,
            ann_layout=ann_layout,
            layer_params = layer_params
        )

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
            if i > 1:
                del(spectra)
                del(coordinates)
            l = 0
            spectra = [spectrum() for _ in xrange(buffer)]
            coordinates = coordinate_arrays(buffer)
            for j in xrange(buffer):
                if j % 100 == 0:
                    print "Reading spectrum #: ", i+j
                if filetype == 1:
                    # If the return value is one, we do a check on coordinates and this 
                    # object is outside the wanted boundaries. Therefore neglect and continue
                    if read_spSpec_fitsio(files[i+j], spectra[j], settings) == 1:
                        continue
                if filetype == 2:
                    if read_spec_fitsio(files[i+j], spectra[j], settings) == 1:
                        continue
                if filetype == 3:
                    if read_speclya_fitsio(files[i+j], spectra[j], settings, resid_corr = resid_corr) == 1:
                        continue
                if filetype == 4:
                    read_mockspec_fitsio(files[i+j], spectra[j])
                if settings.flux_corr == 1:
                    if apply_flux_correction(spectra[j], settings.flux_corr_list[i+j], settings, resid_corr) == 1:
                        spectra[j].z = 999
                        spectra[j].flag = 4
                        continue

            from SDSSmodules.SDSScorrections import fill_coordinate_arrays_from_buffer
            # Fill the l and b arrays with the coordinates of the buffer
            # and get the E(B-V) values from the dustmap
            print 'filling coordinate arrays from buffer...'
            fill_coordinate_arrays_from_buffer(coordinates, spectra, dustmap, 0, buffer)

            ###############################################
            ##### Prepare data for ANN, if we use one #####
            ###############################################
            if settings.predict_from_mlp == True:
                print 'Start reading data for theano'
                if '--unpickle_data' in args:
                    start = i
                    end   = i + buffer
                    path_data    = '/mnt/Daten/Uni/SDSS/data/data_predict/datasets_' + str(start) + '_' + str(end) + '_' + str(1936) + '.dat'
                    saved_data   = open(path_data)
                    datasets     = cPickle.load(saved_data)
                    wave_predict = cPickle.load(saved_data)
                    size         = cPickle.load(saved_data)
                    saved_data.close()
                else:
                    datasets, size, size_out, wave_predict = load_SDSS_predict(None, i, i+buffer, filelist=files, wholespec=True, percentile = 100)
                if filetype == 3:
                    predict_set_x, predict_set_y = datasets[0]
                else:
                    # if not speclya files, only set_x is returned
                    predict_set_x = datasets
                print 'compile theano prediction function'
                predict_from_mlp = theano.function(
                    inputs=[index],
                    outputs=[
                        classifier.y_pred
                        #self.classifier.hiddenLayer[0].output,
                    ],#hiddenLayer[0].output,
                    givens={
                        x: predict_set_x[index:(index+1)]
                    }
                )
#            create_log_arrays(spectra)
        if i % 100 == 0:
            print "Working on spectrum #: ", i
        spectra[l].filename = file
        # Conditions on the QSOs:
        if spectra[l].z >= 2.2 and spectra[l].z <= 5.3:
            # Dust corrections. Only done, if --dust flag is set on startup
            if settings.dust == 1:
                Gal_extinction_correction(spectra[l])

            ###################################
            # Set 'powerlaw' attribute to ANN prediction:
            ###################################
            if settings.predict_from_mlp == True:
                spectra[l].powerlaw = np.zeros(np.size(spectra[l].wave))
                npix_predict        = np.size(wave_predict[l])

                # define the same wavelength limits as used for training (in the end
                # we need a 'powerlaw' array of the same size as the flux array, thus we
                # need to put the prediction into an array of that size
                lya_low    = 1030
                lya_up     = 1585
                # temporary wavelength array in restframe
                wave_temp  = spectra[l].wave / (1.0 + spectra[l].z)
                # get same indices as used for training
                index_wave = np.where((wave_temp > lya_low) &
                                      (wave_temp < lya_up ) &
                                      (np.isfinite(spectra[l].flux) == True))[0]
                size_wave  =  np.size(index_wave)
                # get regulator to get prediction to correct magnitude again
                if index_wave != []:
                    regulator = np.percentile(spectra[l].flux[index_wave], 100)
                else:
                    print 'regulator broken!!'
                    regulator = 1
                # get prediction and apply regulator
                continuum_prediction = predict_from_mlp(l)[0][0] * regulator
                # rebin to same grid as used in spectrum's wavelength array
                continuum_prediction = rebin_1d(continuum_prediction, size_wave, method='spline')
                # put prediction into empty array of size np.size(spec.flux)
                spectra[l].powerlaw[index_wave] = continuum_prediction
                # 'fake' alpha to pass check below:
                spectra[l].alpha = 1
            else:
                ########################################
                # IF we don't use a ANN, fit a powerlaw:
                ########################################
                #powerlaw function
                # fit_powerlaw(spectra[i], 0)
                zem_index = calc_zem_index(spectra[l].z, z_min, z_delta)
                fit_powerlaw_individual(spectra[l], settings, return_data=0, zem_index = zem_index)
                alpha_of_files_used.append(spectra[l].alpha)
                
                #fit_powerlaw_plus_compspec(spectra[i], settings, return_data=0, zem_index = zem_index)
                #spectra[i].powerlaw = spectra[i].model_sdss
                # if spectra[i].PCA_qual == 1:
                #     spectra[i].alpha = 1.0
                # else:
                #    print spectra[i].PCA_qual


            if spectra[l].alpha == -999:
                print 'fitting failed!'
                print i, spectra[l].alpha, spectra[l].z

            #spectra[l].alpha = -999
            # alpha cut
            if spectra[l].alpha < alpha_top and spectra[l].alpha > alpha_low:
                # calculate comp spec
                if build_compspec(compspec, spectra[l]) == 0:
                    compspec.spectra_count += 1
                    # flag spectrum as used and contributed
                    spectra[l].flag = 0
                    files_used.append(file)
                else:
                    files_not_used_compspec.append(file)
                    # flag spectrum as not used, discarded in build_compspec
                    spectra[l].flag = 3
            else:
                files_not_used_alpha.append(file)
                alpha_of_files_not_used.append(spectra[l].alpha)
                # flag spectrum as not used, because of cut on spectral index
                spectra[l].flag = 2
        else:
            files_not_used_zem.append(file)
            # flag spectrum as not used, because of redshift
            if spectra[l].z == 999:
                # means it has bad redshift
                print 'ha?'
                spectra[l].flag == 11
            else:
                spectra[l].flag = 1
        # Every 50th loop, we free all objects, which are not used anymore. 
        # This function is called automatically, but not often enough. Reduces
        # memory usage quite a lot.
        if spectra[l].alpha == -999:
            # print ""
            # print ""
            # print file, spectra[l].alpha, spectra[l].z
            alpha_wrong_count += 1
        if i % 500 == 0:
            gc.collect()
        # Free all big arrays, which won't be needed anymore, after this loop. 
        # Unecessary memory usage.
        del(spectra[l].flux)
        del(spectra[l].flux_error)
        del(spectra[l].wave)
        del(spectra[l].powerlaw)
        del(spectra[l].status)
        del(spectra[l].snr)
        # increase loop variable, since we run over files and not an integer
        l += 1


    ########################################
    ######### cPickle alpha values #########
    ########################################

    alpha_values = open('./alpha_of_files_used.dat', 'wb')
    cPickle.dump(alpha_of_files_used, alpha_values, -1)
    alpha_values.close()

	        
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
