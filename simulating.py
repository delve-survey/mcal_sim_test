import logging
import shutil
import tempfile

import numpy as np
import yaml
import joblib
import galsim
import fitsio
from esutil.ostools import StagedOutFile

from .files import (
    get_band_info_file,
    make_dirs_for_file,
    get_truth_catalog_path,
    expand_path)
from .constants import MEDSCONF
from .truthing import make_coadd_grid_radec
from .sky_bounding import get_rough_sky_bounds, radec_to_uv
from .wcsing import get_esutil_wcs, get_galsim_wcs
from .galsiming import render_sources_for_image
from .psf_wrapper import PSFWrapper

logger = logging.getLogger(__name__)


#reminder for Lucas, Dhayaa and Chihway: lines commented with "CONFIG_change" have hardcoded things that we might want to actually have in a CONFIG file


class End2EndSimulation(object):
    """An end-to-end DES Y3 simulation.
    Parameters
    ----------
    seed : int
        The seed for the global RNG.
    output_meds_dir : str
        The output DEADATA/MEDS_DIR for the simulation data products.
    tilename : str
        The DES coadd tile to simulate.
    bands : str
        The bands to simulate.
    gal_kws : dict
        Keyword arguments to control the galaxy content of the simulation.
        Right now these should include:
            n_grid : int
                The galaxies will be put on a grid with `n_grid`
                on a side.
            g1 : float
                The true shear on the one-axis.
            g2 : float
                The true shear on the two-axis.
    psf_kws : dict
        Kyword arguments to control the PSF used for the simulation.
        Right now these should include:
            type : str
                One of 'gauss' and that's it.
    Methods
    -------
    run()
        Run the simulation, writing the data to disk.
    """
    def __init__(self, *,
                 seed, output_meds_dir, tilename, bands,
                 gal_kws, psf_kws):
        self.output_meds_dir = output_meds_dir
        self.tilename = tilename
        self.bands = bands
        self.gal_kws = gal_kws
        self.psf_kws = psf_kws
        self.seed = seed
        # any object within a 128 coadd pixel buffer of the edge of a CCD
        # will be rendered for that CCD
        self.bounds_buffer_uv = 128 * 0.263

        self.draw_method = 'auto'

        # make the RNGS
        # one for galaxies in the truth catalog
        # one for noise in the images
        seeds = np.random.RandomState(seed=seed).randint(
            low=1, high=2**30, size=2)
        self.gal_rng = np.random.RandomState(seed=seeds[0])
        self.noise_rng = np.random.RandomState(seed=seeds[1])

        # load the image info for each band
        self.info = {}
        for band in bands:
            fname = get_band_info_file( #this function is picking out the filename
                meds_dir=self.output_meds_dir,
                medsconf=MEDSCONF,
                tilename=self.tilename,
                band=band)
            with open(fname, 'r') as fp:
                self.info[band] = yaml.load(fp, Loader=yaml.Loader)

    def run(self):
        """Run the simulation w/ galsim, writing the data to disk."""

        logger.info(' simulating coadd tile %s', self.tilename)

        # step 1 - make the truth catalog
        truth_cat = self._make_truth_catalog() #"TRUTH" because it contains the input positions of future simulated galaxies - as opposed to sextracting galaxies

        # step 2 - per band, write the images to a tile
        for band in self.bands:
            self._run_band(band=band, truth_cat=truth_cat)

    def _run_band(self, *, band, truth_cat):
        """Run a simulation of a truth cat for a given band."""

        logger.info(" rendering images in band %s", band)

        noise_seeds = self.noise_rng.randint(
            low=1, high=2**30, size=len(self.info[band]['src_info']))

        jobs = []
        for noise_seed, se_info in zip(
                noise_seeds, self.info[band]['src_info']):

            src_func = LazySourceCat(
                truth_cat=truth_cat,
                wcs=get_galsim_wcs(
                    image_path=se_info['image_path'],
                    image_ext=se_info['image_ext']),
                psf=self._make_psf_wrapper(se_info=se_info), #STOPPED HERE, so continue from _make_PSF_wrapper
                g1=self.gal_kws['g1'],
                g2=self.gal_kws['g2'])

            jobs.append(joblib.delayed(_render_se_image)(
                se_info=se_info,
                band=band,
                truth_cat=truth_cat,
                bounds_buffer_uv=self.bounds_buffer_uv,
                draw_method=self.draw_method,
                noise_seed=noise_seed,
                output_meds_dir=self.output_meds_dir,
                src_func=src_func))

        with joblib.Parallel(
                n_jobs=-1, backend='loky', verbose=50, max_nbytes=None) as p:
            p(jobs)

    def _make_psf_wrapper(self, *, se_info):
        if self.psf_kws['type'] == 'gauss':
            psf_model = galsim.Gaussian(fwhm=0.9)
        elif self.psf_kws['type'] == 'piff':
            from ..des_piff import DES_Piff
            psf_model = DES_Piff(expand_path(se_info['piff_path']))
            assert self.draw_method == 'auto'
        elif self.psf_kws['type'] == 'gauss-pix':
            from .gauss_pix_psf import GaussPixPSF
            kwargs = {k: self.psf_kws[k] for k in self.psf_kws if k != 'type'}
            psf_model = GaussPixPSF(**kwargs)
            assert self.draw_method == 'auto'
        else:
            raise ValueError(
                "psf type '%s' not recognized!" % self.psf_kws['type'])

        wcs = get_galsim_wcs(
            image_path=se_info['image_path'],
            image_ext=se_info['image_ext'])
        psf_wrap = PSFWrapper(psf_model, wcs)

        return psf_wrap

    def _make_truth_catalog(self):
        """Make the truth catalog."""
        # always done with first band
        band = self.bands[0]
        coadd_wcs = get_esutil_wcs( #gets WCS and seems to cache it into memory
            image_path=self.info[band]['image_path'],
            image_ext=self.info[band]['image_ext'])

        ra, dec, x, y = make_coadd_grid_radec( #gets source positions on a grid, presumably over a sky footprint and dithers them, for eg. simulations without blending
            rng=self.gal_rng, coadd_wcs=coadd_wcs,
            return_xy=True, n_grid=self.gal_kws['n_grid'])

        truth_cat = np.zeros(
            len(ra), dtype=[
                ('number', 'i8'),
                ('ra', 'f8'),
                ('dec', 'f8'),
                ('x', 'f8'),
                ('y', 'f8')])
        truth_cat['number'] = np.arange(len(ra)).astype(np.int64) + 1
        truth_cat['ra'] = ra
        truth_cat['dec'] = dec
        truth_cat['x'] = x
        truth_cat['y'] = y

        truth_cat_path = get_truth_catalog_path( #literally just joins strings
            meds_dir=self.output_meds_dir,
            medsconf=MEDSCONF,
            tilename=self.tilename)

        make_dirs_for_file(truth_cat_path)
        fitsio.write(truth_cat_path, truth_cat, clobber=True)

        return truth_cat

##############################
# from https://github.com/beckermr/misc/blob/main/matts_misc/simple_des_y3_sims/files.py

def get_band_info_file(*, meds_dir, medsconf, tilename, band):
    """Get the path of the YAML file holding the info dict for the
    `tilename` and `band`.
    Parameters
    ----------
    meds_dir : str
        The DESDATA/MEDS_DIR path where the info file is located.
    medsconf : str
        The MEDS file version (e.g., 'y3v02').
    tilename : str
        The DES coadd tilename (e.g., 'DES2122+0001').
    bands : str
        A bands (e.g., 'r').
    Returns
    -------
    info_file : str
        The YAML file with the coadd + SE information.
    """
    return os.path.join(
        meds_dir,
        'simple_des_y3_sims',
        medsconf,
        'band_info_files',
        '%s_%s_info.yaml' % (tilename, band))


@lru_cache(maxsize=256)
def get_esutil_wcs(*, image_path, image_ext):
    """Read the WCS solution of an image into an esutil.wcsutil.WCS object.
    Parameters
    ----------
    image_path : str
        The path to the image.
    image_ext : int or str
        The extension with the WCS information.
    Returns
    -------
    wcs : esutil.wcsutil.WCS
        The WCS object.
    """
    hd = fitsio.read_header(
        image_path, ext=image_ext)
    hd = {k.lower(): hd[k] for k in hd if k is not None}
    return esutil.wcsutil.WCS(hd)

def make_coadd_grid_radec(*, n_grid, coadd_wcs, rng, return_xy=False): #creates a grid of simulations, positions in the grid are dithered here
    """Make a grid of points in the coadd image coordinate system and
    return their locations in ra-dec.
    Parameters
    ----------
    n_grid : int
        The number of objects across the grid in each direction. The total
        number of objects will be `n_grid**2`.
    coadd_wcs : esutil.wcsutil.WCS
        The coadd WCS solution.
    rng : np.random.RandomState
        An RNG to use. This RNg is used to dither the locations on the coadd
        grid within a pixel.
    return_xy : bool, optional
        If True, also return the x and y positions. Default is False
    Returns
    -------
    ra : np.ndarray
        The array of ra positions of the sources.
    dec : np.ndarray
        The array of dec positions of the sources.
    x : np.ndarray
        The array of column positions. Only returned if `return_xy=True`.
    y : np.ndarray
        The array of row positions. Only returned if `return_xy=True`.
    """
    L = 10000  # hard code this since it will not change
    dL = L / n_grid
    dL_2 = dL / 2

    x = []
    y = []
    for row_ind in range(n_grid):
        for col_ind in range(n_grid):
            _x = col_ind * dL + dL_2 + 1
            _y = row_ind * dL + dL_2 + 1

            # dither
            _x += rng.uniform(low=-0.5, high=0.5)
            _y += rng.uniform(low=-0.5, high=0.5)

            x.append(_x)
            y.append(_y)

    x = np.array(x)
    y = np.array(y)
    ra, dec = coadd_wcs.image2sky(x, y)

    if return_xy:
        return ra, dec, x, y
    else:
        return ra, dec


def get_truth_catalog_path(*, meds_dir, medsconf, tilename):
    """Get the truth catalog path.
    Parameters
    ----------
    meds_dir : str
        The DESDATA/MEDS_DIR path where the info file is located.
    medsconf : str
        The MEDS file version (e.g., 'y3v02').
    tilename : str
        The DES coadd tilename (e.g., 'DES2122+0001').
    Returns
    -------
    truth_file_path : str
        The path to the truth file.
    """
    return os.path.join(
        meds_dir,
        'simple_des_y3_sims',
        medsconf,
        'truthcats',
        '%s_truthcat.fits' % tilename
    )

def make_dirs_for_file(filename):
    """Make all of the parent directories for a file at `filename`."""
    dirname = os.path.dirname(filename)
    if len(dirname) > 0:
        os.makedirs(dirname, exist_ok=True)

class LazySourceCat(object):
    """A lazy source catalog that only builds objects to be rendered as they
    are needed.
    Parameters
    ----------
    truth_cat : structured np.array
        The truth catalog as a structured numpy array.
    wcs : galsim.GSFitsWCS
        A galsim WCS instance for the image to be rendered.
    psf : PSFWrapper
        A PSF wrapper object to use for the PSF.
    g1 : float
        The shear to apply on the 1-axis.
    g2 : float
        The shear to apply on the 2-axis.
    Methods
    -------
    __call__(ind)
        Returns the object to be rendered from the truth catalog at
        index `ind`.
    """
    def __init__(self, *, truth_cat, wcs, psf, g1, g2):
        self.truth_cat = truth_cat
        self.wcs = wcs
        self.psf = psf
        self.g1 = g1
        self.g2 = g2

    def __call__(self, ind):
        pos = self.wcs.toImage(galsim.CelestialCoord(
            ra=self.truth_cat['ra'][ind] * galsim.degrees,
            dec=self.truth_cat['dec'][ind] * galsim.degrees))
        obj = galsim.Exponential(half_light_radius=0.5).withFlux(64000) #might replace these hardcoded numbers by some CONFIG_change file input 
        obj = obj.shear(g1=self.g1, g2=self.g2)
        psf = self.psf.getPSF(image_pos=pos)
        return galsim.Convolve([obj, psf]), pos

@lru_cache(maxsize=256)
def get_galsim_wcs(*, image_path, image_ext):
    """Read the WCS solution of an image into a galsim WCS object.
    Parameters
    ----------
    image_path : str
        The path to the image.
    image_ext : int or str
        The extension with the WCS information.
    Returns
    -------
    wcs : galsim WCS
        The WCS object.
    """
    hd = fitsio.read_header(
        image_path, ext=image_ext)
    hd = {k.upper(): hd[k] for k in hd if k is not None}
    wcs = galsim.FitsWCS(header=hd)
    assert not isinstance(wcs, galsim.PixelScale)  # this has been a problem
    return wcs
