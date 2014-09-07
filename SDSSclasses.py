# contains astropy coordinate objects used to work with coordinates
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy as hp
import numpy as np


################################################################################
############################## Classes #########################################
################################################################################

# TODO: eventually change class names to capital letters!
# use @property for arrays and have one big numpy array for flux, flux_error etc?
# would need to move creation of arrays from main function into read function though,
# since need to provide length of arrays

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
        self.alpha       = -999         # spectral index of spectrum
        self.alpha_error = -999         # 1 sigma error of alpha
        self.beta        = -999         # alpha = -beta - 2
        self.delta       = -999         # start of powerlaw wavelength array. 
                                        # bad name, but taken from C code
        self.chisq       = -999         # reduced chi^2 from fit to log log data
        # Coordinates
        # using astropy.coordinates object. Currently using v0.3.2 of astropy
        # from v0.4 (in develeopment atm) on exchanged with SkyCoord. 
        # Done.
        self.coordinates = SkyCoord(ra = 0*u.degree, dec = 0*u.degree, frame='fk5')

        self.npix           = 0         # number of pixels in spectrum
        # Arrays:
        self.flux           = []        # flux of spectrum
        self.flux_error     = []        # 1 sigma of flux
        self.powerlaw       = []        # powerlaw continuum
        self.powerlaw_error = []        # 1 sigma of p.c.
        self.model_sdss     = []        # PCA model by SDSS collaboration (continuum)
        self.flux_norm      = []
        self.flux_norm_error= []
        self.wave           = []        # wavelength 
        self.emfree         = []        # Emission free regions used for fitting
        self.mag            = []        # magnitudes in ugriz band
        self.mag_error      = []        # 1 sigma of magnitudes
        self.smag           = []        # spectroscopic magnitudes
        self.status         = []        # status; 1 = good pixel, 0 = bad pixel
        self.snr            = []        # signal to noise ratio
        self.flag           = 0         # flag, which flags if the quasar was:
                                        # 0  =  used, contributed to cspec
                                        # 1  =  not used, because failed red shift check
                                        # 2  =  not used, because failed spectral index cut
                                        # 3  =  not used, discarded in build_compspec

        self.filename = []              # Filename of FITS file QSO is stored in
        self.plateid = 0                # PlateID of spectrum
        self.fiberid = 0                # FiberID of spectrum
        self.MJD = 0                    # MJD of spectrum taken on

# Class for the 'composite spectrum'
# this object contains the opacity in the flux array after SDSScompspec was used
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
        self.powerlaw      = []         # powerlaw continuum fit to sum
        self.normalized    = []         # normalized composite spectrum
        self.spectra_count = 0          # number of used spectras in total
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
            self.wave.append(10.0**(3.57520 + 0.0001*i))
#            self.wave.append(10.0**(3.5899 + 0.0001*i))
#            self.wave.append(10.0**(3.58020 + 0.0001*i))

# Object for real composite spectrum
class composite_spectrum:
    def __init__(self, npix, i):
        self.wave          = []         # wavelength
        self.flux          = []         # flux of composite spectrum
        self.flux_error    = []         # 1 sigma of opacity
        self.sum           = []         # flux_i / powerlaw_i sum over all i QSOs 
        self.sum2          = []         # sum of squares of the above
        self.nhist         = []         # number of contributing spectra to each
                                        # element in sum and sum2
        self.powerlaw      = []         # powerlaw continuum fit to sum
        self.normalized    = []         # normalized composite spectrum
        self.npix          = npix
        self.spectra_count = 0          # number of used spectras for this compspec
        self.zem_index     = i          # index of redshift bin
        self.min_z_in_bin  = 0          # lowest emission redshift of object contributed
        self.max_z_in_bin  = 0          # highest emission redshift of object contributed
        self.zem_in_bin    = []         # emission redshifts contributed
        self.mean_a        = 0          # mean of alpha of used spectra
        self.sigma_a       = 0          # 1 sigma of alpha of used spectra
        self.median_a      = 0          # median of alpha of used spectra
        self.siqr_a        = 0          # semi-interquartile range of alpha
        for i in xrange(npix):
            # initialize sum, sum2 and nhist to 0 and wavelength array starting
            # from 900 Angstrom
            self.sum.append(0)          
            self.sum2.append(0)
            self.nhist.append(0)
            self.wave = 10**(2.9542 + 0.0001*np.arange(npix))


# Class to store coordinate arrays
class coordinate_arrays:
    def __init__(self, nspec):
        self.ra_array  = np.zeros(nspec)
        self.dec_array = np.zeros(nspec)
        self.l_array   = np.zeros(nspec)
        self.b_array   = np.zeros(nspec)

# Class for the program settings

class program_settings:
    def __init__(self):
        self.dust             = 0
        self.nprocs           = 8
        self.program_name     = 'PyS_SDSScompspec'
        self.inputfile        = ''
        self.outfile          = ''
        self.outfile_basis    = ''
        self.spectra_list     = ''
        self.cspec            = 0
        self.deviation_factor = 3.0
        self.z_max            = 8.0  # maximum emission redshift for creation
                                     # of redshift bins
        self.z_min            = 2.0  # minimum emission redshift for creation
                                     # of redshift bins
        self.z_delta          = 0.4  # delta /width of redshift bins 
        self.z_bins           = int((self.z_max - self.z_min)/self.z_delta)
        self.coords           = 0
        # Using HEALpy, we have no use for the coordinate specific variables
        # self.l_min        = -999
        # self.l_max        = -999
        # self.b_min        = -999
        # self.b_max        = -999
        # self.delta_l      = 10
        # self.delta_b      = self.delta_l / 2.0
        # Instead, we determine the resolution of the map by the HEALpix parameter
        # Nside. Default: 8, for runtime while testing)
        self.map_nside        = 2
        self.map_npix         = hp.nside2npix(self.map_nside)
