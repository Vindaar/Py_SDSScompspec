# contains astropy coordinate objects used to work with coordinates
from astropy.coordinates import ICRS, Galactic
from astropy import units as u


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
# TODO: Think about switching to numpy arrays
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
        
