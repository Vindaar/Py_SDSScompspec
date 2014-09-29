# This module contains the functions, which deal with the composite spectra
# objects. 




################################################################################
############################## Compspec ########################################
################################################################################
    
def build_compspec(cspec, spec):
# this function calculates the 'composite spectrum' from the individual spectra
# by adding flux/powerlaw values in a certain Restframe wavelength range
    # Define Restrange array:
    Restrange = []
    # 1095, 1150
    # 1110, 1135 tested
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
    iterator = int(find_element_larger_in_arrays(spec.wave, cspec.wave[0], spec.npix))
    if spec.wave[iterator] > cspec.wave[0] and iterator > 0:
        iterator -= 1
    citerator = int(find_element_larger_in_arrays(cspec.wave, spec.wave[0], spec.npix))

    params_iter = open("params_iter.txt", "a")
    params_iter.write(str(iterator)+"    "+str(citerator)+"     "+str(spec.filename))

    # experimental:
    # Does not give exact results as old version below, but MUCH faster!
    index = np.where((spec.wave                 != 0)            &
                     (spec.wave                 > cspec.wave[0]) &
                     (spec.wave                 > spec.lambmin)  &
                     (spec.wave                 < spec.lambmax-10**(0.0001))  &
                     (spec.mask_comb            == 0)            &
                     (np.isnan(spec.flux)       == False)        &
                     (np.isnan(spec.flux_error) == False)        &
                     (spec.flux_error           != 0)            &
                     (spec.powerlaw             != 0))[0]

    indices_for_cspec = np.subtract(index, iterator)
    indices_for_cspec = np.add(indices_for_cspec, citerator)#np.arange(citerator, spec.npix - iterator)
    indices_for_spec  = index#np.add(index, iterator)

    temp_sum   = np.zeros(5763)
    temp_sum2  = np.zeros(5763)
    temp_nhist = np.zeros(5763)
    np.put(temp_sum,   indices_for_cspec, spec.flux[indices_for_spec] / spec.powerlaw[indices_for_spec])
    np.put(temp_sum2,  indices_for_cspec, (spec.flux[indices_for_spec] / spec.powerlaw[indices_for_spec])**2)
    np.put(temp_nhist, indices_for_cspec, 1)

#    print temp_sum
#    print spec.flux[index] / spec.powerlaw[index]

    cspec.sum   = np.add(cspec.sum,   temp_sum)
    cspec.sum2  = np.add(cspec.sum2,  temp_sum2)
    cspec.nhist = np.add(cspec.nhist, temp_nhist)


    # MUCH slower (50ms vs. 500mu s per call) than above. However, slightly different results
    # Could compile in Cython and see how fast then
    # Run over all pixels / wavelengths of the spectrum
    # - iterator, because we don't want to access elements outside of array bounds
    # from 1 to pixels - iterator - 1, because we check for lambmin and lambmax
    # by looking at element i - 1, i + 1 respectively
    # for i in xrange(1, spec.npix - iterator - 1):
    #     # if a lot of checks are true, 
    #     if(spec.wave[i+iterator]                             != 0             and
    #        spec.wave[i+iterator] + spec.wave[i+iterator - 1] > 2*spec.lambmin and
    #        spec.wave[i+iterator] + spec.wave[i+iterator + 1] < 2*spec.lambmax and
    #        np.isnan(spec.flux[i+iterator])                   == 0             and
    #        np.isnan(spec.flux_error[i+iterator])             == 0             and
    #        spec.flux_error[i+iterator]                       != 0             and
    #        spec.powerlaw[i+iterator]                         != 0):
    #         # we actually use the spectrum to add to the composite spectrum
    #         cspec.sum[i+citerator]   += spec.flux[i+iterator]/spec.powerlaw[i+iterator]
    #         cspec.sum2[i+citerator]  += (spec.flux[i+iterator]/spec.powerlaw[i+iterator])**2
    #         # and we count this contribution to this wavelength
    #         cspec.nhist[i+citerator] += 1
    #     else: 
    #         # if checks not met, go to next pixel
    #         continue
    #     # if we would be outside of lambmax in the next iteration, stop for loop
    #     if spec.wave[i+iterator] + spec.wave[i+iterator + 1] > 2*spec.lambmax: 
    #         break

    params_iter.close()
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

    parameters = open('parameters.txt', 'a')
    # eventually take parameters out again. No real use except debugging

    for i in xrange(np.size(cspec.sum)):
        parameters.write(str(cspec.sum[i])+"     "+str(cspec.nhist[i])+"     "+str(cspec.sum2[i])+'\n')
    parameters.close()

    # Changed 3000 pixels to 4000, because in DR10 big sidx values appear
    for i in xrange(4000):
        # if nhist is non zero (because we divide by it) 
        if cspec.nhist[i] > 0:
            # we add the flux of the composite spectrum as:
            # sum / nhist
            # TODO: for some reason some values of the flux (around i = 2700) are below 0. Shouldn't be
            # the case I think?
            # TODO: Maybe change to numpy arrays? 
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

  
