# This module contains the color correction for the PSF magnitudes
# we don't use this module currently


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
    wavelength = 10.0**(beginwl + 0.0001*np.arange(npix))
#    for i in xrange(npix):
#        wavelength.append(10.0**(beginwl + 0.0001*i))
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

    # Calculate the sum of the color bands fluxes
    # Commented code below is more obvious, but slower
    usum = sum(a[j:spec.npix+j,0]*spec.flux)
    gsum = sum(a[j:spec.npix+j,1]*spec.flux)
    rsum = sum(a[j:spec.npix+j,2]*spec.flux)
    isum = sum(a[j:spec.npix+j,3]*spec.flux)
    zsum = sum(a[j:spec.npix+j,4]*spec.flux)

    # while i < spec.npix and j < npix:
    #     usum += a[j][0] * spec.flux[i]
    #     gsum += a[j][1] * spec.flux[i]
    #     rsum += a[j][2] * spec.flux[i]
    #     isum += a[j][3] * spec.flux[i]
    #     zsum += a[j][4] * spec.flux[i]
    #     i += 1
    #     j += 1
        
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
