# contains astropy coordinate objects used to work with coordinates
from astropy.coordinates import ICRS, Galactic
from astropy import units as u


################################################################################
############################## Classes #########################################
################################################################################
 
# Class for the spectra
class spectrum:
    def __init__(self):
        self.beginwl     = 0            # wavelength of first pixel in spectrum
        self.deltawl     = 0            # increment of wavelength per pixel 
        # lower and upper bound of interval used for compspec scaled to QSOs 
        # emission redshift: (1+z)*[1095;1150]
        self.lambmin     = 0            
        self.lambmax     = 0
        self.cpix        = 0            # 1
        self.z           = 0            # emission redshift of the QSO
        # Fitting parameters
        self.alpha       = 0            # spectral index of spectrum
        self.alpha_error = 0            # 1 sigma error of alpha
        self.beta        = 0            # alpha = -beta - 2
        self.delta       = 0            # start of powerlaw wavelength array. 
                                        # bad name, but taken from C code
        # Coordinates
        # using astropy.coordinates object. Currently using v0.3.2 of astropy
        # from v0.4 (in develeopment atm) on exchanged with SkyCoord. 
        self.coordinates = ICRS(ra = 0, dec = 0, unit=(u.degree, u.degree))

        self.npix           = 0         # number of pixels in spectrum
        # Arrays:
        self.flux           = []        # flux of spectrum
        self.flux_error     = []        # 1 sigma of flux
        self.powerlaw       = []        # powerlaw continuum
        self.powerlaw_error = []        # 1 sigma of p.c.
        self.flux_norm      = []
        self.flux_norm_error= []
        self.wave           = []        # wavelength 
        self.emfree         = []        # Emission free regions used for fitting
        self.mag            = []        # magnitudes in ugriz band
        self.mag_error      = []        # 1 sigma of magnitudes
        self.smag           = []        # spectroscopic magnitudes
        self.status         = []        # status; 1 = good pixel, 0 = bad pixel
        self.snr            = []        # signal to noise ratio

        self.filename = []              # Filename of FITS file QSO is stored in
        self.plateid = 0                # PlateID of spectrum
        self.fiberid = 0                # FiberID of spectrum
        self.MJD = 0                    # MJD of spectrum taken on

# Class for the composite spectrum
# TODO: Think about switching to numpy arrays
class comp_spectrum:
    def __init__(self, npix):
        self.wave          = []         # wavelength
        self.flux          = []         # actually opacity not flux
        self.flux_error    = []         # 1 sigma of opacity
        self.sum           = []         # flux_i / powerlaw_i sum over all i QSOs 
        self.sum2          = []         # sum of squares of the above
        self.nhist         = []         # number of contributing spectra to each
                                        # element in sum and sum2
        self.spectra_count = 0          # number of used specrtas in total
        self.mean_a        = 0          # mean of alpha of used spectra
        self.sigma_a       = 0          # 1 sigma of alpha of used spectra
        self.median_a      = 0          # median of alpha of used spectra
        self.siqr_a        = 0          # semi-interquartile range of alpha
        for i in xrange(npix):
            # initialize sum, sum2 and nhist to 0 and wavelength array starting
            # from 10**3.58020 Angstrom
            self.sum.append(0)          
            self.sum2.append(0)
            self.nhist.append(0)
            self.wave.append(10.0**(3.58020 + 0.0001*i))

# Class for the program settings

class program_settings:
    def __init__(self):
        self.dust         = 0
        self.nprocs       = 8
        self.program_name = 'PyS_SDSScompspec'
        self.inputfile    = ''
        self.outfile      = ''
        self.cspec        = 0
        self.coords       = 0
        self.l_min        = -999
        self.l_max        = -999
        self.b_min        = -999
        self.b_max        = -999
#        self.coords_min = ICRS(ra = 0, dec = 0, unit=(u.degree, u.degree))
#        self.coords_max = ICRS(ra = 0, dec = 0, unit=(u.degree, u.degree))
        
