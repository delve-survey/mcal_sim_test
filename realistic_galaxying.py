import os
import functools
import collections
import numpy as np

import galsim
import fitsio

WLDeblendData = collections.namedtuple(
    'WLDeblendData',
    [
        'cat', 'rand_rot', 'survey_name', 'bands', 'surveys',
        'builders', 'total_sky', 'noise', 'ngal_per_arcmin2',
        'psf_fwhm', 'pixel_scale',
    ],
)

CosmosData = collections.namedtuple(
    'CosmosData',
    [
        'cat', 'rand_rot',
    ],
)


# @functools.lru_cache(maxsize=8)
def _cached_catalog_read():
    fname = os.path.join(os.environ.get('CATSIM_DIR', '.'), 'OneDegSq.fits',)
    return fitsio.read(fname)


def _cached_COSMOS_catalog_read():
    fname = os.path.join(os.environ.get('CATCOSMOS_DIR', '.'), 'input_cosmos_v4.fits',)
    return fitsio.read(fname)


# @functools.lru_cache(maxsize=8)
def init_descwl_catalog(*, survey_bands, rng):
    """Initialize weak lensing deblending survey data.

    Parameters
    ----------
    survey_bands : str
        The name of the survey followed by the bands like 'des-riz', 'lsst-iz', etc.
        
    rng : np.random.default_rng instance
        The rng instance that will be used to generate the random rotation
        of the galaxy

    Returns
    -------
    data : WLDeblendData
        Namedtuple with data for making galaxies via the weak lesning
        deblending package.
    """
    survey_name, bands = survey_bands.split("-")
    bands = [b for b in bands]

    if survey_name not in ["des", "lsst"]:
        raise RuntimeError(
            "Survey for wldeblend must be one of 'des' or 'lsst'"
            " - got %s!" % survey_name
        )

    if survey_name == "lsst":
        scale = 0.2
    elif survey_name == "des":
        scale = 0.263

    # guard the import here
    import descwl

    # set the exposure times
    if survey_name == 'des':
        exptime = 90 * 10
    else:
        exptime = None

    wldeblend_cat = _cached_catalog_read()

    surveys = []
    builders = []
    total_sky = 0.0
    for iband, band in enumerate(bands):
        # make the survey and code to build galaxies from it
        pars = descwl.survey.Survey.get_defaults(
            survey_name=survey_name.upper(),
            filter_band=band)

        pars['survey_name'] = survey_name
        pars['filter_band'] = band
        pars['pixel_scale'] = scale

        # note in the way we call the descwl package, the image width
        # and height is not actually used
        pars['image_width'] = 100
        pars['image_height'] = 100

        # reset the exposure times if we want
        if exptime is not None:
            pars['exposure_time'] = exptime

        # some versions take in the PSF and will complain if it is not
        # given
        try:
            _svy = descwl.survey.Survey(**pars)
        except Exception:
            pars['psf_model'] = None
            _svy = descwl.survey.Survey(**pars)

        surveys.append(_svy)
        builders.append(descwl.model.GalaxyBuilder(
            survey=surveys[iband],
            no_disk=False,
            no_bulge=False,
            no_agn=False,
            verbose_model=False))

        total_sky += surveys[iband].mean_sky_level

    noise = np.sqrt(total_sky)

    if survey_name == "lsst":
        psf_fwhm = 0.85
    elif survey_name == "des":
        psf_fwhm = 1.1

    # when we sample from the catalog, we need to pull the right number
    # of objects. Since the default catalog is one square degree
    # and we fill a fraction of the image, we need to set the
    # base source density `ngal`. This is in units of number per
    # square arcminute.
    ngal_per_arcmin2 = wldeblend_cat.size / (60 * 60)

    
    #CUT OUT LARGE GALAXIES FROM DATASET
    #Check largest axis size and remove galaxy based on that size
#     size = np.max([wldeblend_cat['a_d'], wldeblend_cat['b_d'], wldeblend_cat['a_b'], wldeblend_cat['b_b']], axis = 0)
    size = np.max([wldeblend_cat['a_d'],wldeblend_cat['a_b']], axis = 0)
#     wldeblend_cat['size'] = size
#     wldeblend_cat = wldeblend_cat[size < 0.8]
    wldeblend_cat = wldeblend_cat[size < 0.5]
    
    #If rng not supplied then don't do random rotation
    if rng is None:
        angle = None
    else:
        angle = rng.uniform(low = 0, high = 1, size = len(wldeblend_cat))*360
    
    return WLDeblendData(
        wldeblend_cat, angle, survey_name, bands, surveys,
        builders, total_sky, noise, ngal_per_arcmin2,
        psf_fwhm, scale,
    )


def init_cosmos_catalog(*, rng):
    
    
    cosmos_cat = _cached_COSMOS_catalog_read()
    
    #If rng not supplied then don't do random rotation
    if rng is None:
        angle = None
    else:
        angle = rng.uniform(low = 0, high = 1, size = len(cosmos_cat))*360
        
    return CosmosData(cosmos_cat, angle)


def get_descwl_galaxy(*, descwl_ind, rng, data):
    """Draw a galaxy from the weak lensing deblending package.

    Parameters
    ----------
    descwl_ind : int
        Index of galaxy in descwl catalog. Needed so galaxy in
        every band/exposure looks the same.
    rng : np.random.RandomState
        An RNG to use for galaxy orientation
    data : WLDeblendData
        Namedtuple with data for making galaxies via the weak lesning
        deblending package.

    Returns
    -------
    gal : galsim Object
        The galaxy as a galsim object.
    """
    
    #Overriding rng because we need rotation to be
    #the same for every band for a given galaxy
    
#     rng      = np.random.default_rng(seed = descwl_ind)

#     angle    = rng.uniform() * 360
#     pa_angle = rng.uniform() * 360

#     angle = 0
#     print(data.cat['pa_disk'].flags)
#     print(data.cat['pa_bulge'].flags)
#     data.cat['pa_disk'][rind] = pa_angle
#     data.cat['pa_bulge'][rind] = pa_angle


    return galsim.Sum([
        data.builders[band].from_catalog(
            data.cat[descwl_ind], 0, 0,
            data.surveys[band].filter_band).model.rotate(
                data.rand_rot[descwl_ind] * galsim.degrees)
        for band in range(len(data.builders))
    ])


def get_cosmos_galaxy(*, cosmos_ind, rng, data, band = None):
    """Draw a galaxy from the DES_COSMOS model.

    Parameters
    ----------
    cosmos_ind : int
        Index of galaxy in cosmos catalog. Needed so galaxy in
        every band/exposure looks the same.
    rng : np.random.RandomState
        An RNG to use for galaxy orientation
    data : Cosmos data fitsio file
        Namedtuple with data for making galaxies via the weak lesning
        deblending package.
    band : character
        A single character containing the band of the galaxy.
        If None then average over riz bands.

    Returns
    -------
    gal : galsim Object
        The galaxy as a galsim object.
    """
    
    bulge_frac = data.cat['bdf_fracdev'][cosmos_ind] #Fraction of bulge to total
    
    if band == None:
        flux = np.sum([data.cat['flux_%s'%i][cosmos_ind] for i in ['r', 'i', 'z']])
        print("Why am I in here", band)
    else:
        flux = data.cat['flux_%s'%band][cosmos_ind]
        
    disk  = galsim.Exponential(flux = flux,   half_light_radius = data.cat['bdf_hlr'][cosmos_ind])
    bulge = galsim.DeVaucouleurs(flux = flux, half_light_radius = data.cat['bdf_hlr'][cosmos_ind])
    
    prof  = bulge_frac*bulge + (1 - bulge_frac)*disk
    prof  = prof.shear(g1 = data.cat['bdf_g1'][cosmos_ind], g2 = data.cat['bdf_g2'][cosmos_ind])
    prof  = prof.rotate(data.rand_rot[cosmos_ind]*galsim.degrees)

    return prof


def get_psf_config_wldeblend(*, data):
    """Get a config dict for a the PSF model for the weak lensing deblending
    objects.

    Parameters
    ----------
    data : WLDeblendData
        Namedtuple with data for making galaxies via the weak lesning
        deblending package.

    Returns
    -------
    gs_config : dict
        A dictionary with the PSF info.
    """
    gs_config = {}
    gs_config["type"] = "Kolmogorov"
    gs_config["fwhm"] = data.psf_fwhm
    return gs_config