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
            fname = get_band_info_file(
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
        truth_cat = self._make_truth_catalog()

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
                psf=self._make_psf_wrapper(se_info=se_info),
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
        coadd_wcs = get_esutil_wcs(
            image_path=self.info[band]['image_path'],
            image_ext=self.info[band]['image_ext'])

        ra, dec, x, y = make_coadd_grid_radec(
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

        truth_cat_path = get_truth_catalog_path(
            meds_dir=self.output_meds_dir,
            medsconf=MEDSCONF,
            tilename=self.tilename)

        make_dirs_for_file(truth_cat_path)
        fitsio.write(truth_cat_path, truth_cat, clobber=True)

        return truth_cat