import logging
import numpy as np

import joblib
import esutil as eu
import fitsio
from ngmix import ObsList, MultiBandObsList
from ngmix.gexceptions import GMixRangeError
#from ngmix_compat import NGMixMEDS, MultiBandNGMixMEDS, NGMIX_V1
from ngmix.medsreaders import NGMixMEDS, MultiBandNGMixMEDS#, NGMIX_V1
NGMIX_V1=False #wild guess here

from files import get_meds_file_path, get_mcal_file_path, make_dirs_for_file
from metacal.metacal_fitter import MetacalFitter
from constants import MEDSCONF, MAGZP_REF
from interpolate import interpolate_image_at_mask
import galsim
from run_utils import _run_mcal_one_chunk

logger = logging.getLogger(__name__)

def run_metacal(*, tilename, output_meds_dir, bands, seed, mcal_config):
    """Run metacal on a tile.

    Parameters
    ----------
    tilename : str
        The DES coadd tile on which to run metacal.
    output_meds_dir : str
        The output DEADATA/MEDS_DIR for the simulation data products.
    bands : str
        The bands on which to run metacal.
    seed : int
        The seed for the global RNG.
    mcal_config : yaml file
        The config file for the metacal run
    """
    meds_files = [
        get_meds_file_path(
            meds_dir=output_meds_dir,
            medsconf=MEDSCONF,
            tilename=tilename,
            band=band)
        for band in bands]
    with NGMixMEDS(meds_files[0]) as m:
        cat = m.get_cat()
    logger.info(' meds files %s', meds_files)

    n_chunks = joblib.externals.loky.cpu_count()
    n_obj_per_chunk = cat.size // n_chunks
    if n_obj_per_chunk * n_chunks < cat.size:
        n_obj_per_chunk += 1
    assert n_obj_per_chunk * n_chunks >= cat.size
    logger.info(' running metacal for %d objects in %d chunks', cat.size, n_chunks)

    seeds = np.random.RandomState(seed=seed).randint(1, 2**30, size=n_chunks)

    
    jobs = []
    for chunk in range(n_chunks):
        start = chunk * n_obj_per_chunk
        end = min(start + n_obj_per_chunk, cat.size)
        jobs.append(joblib.delayed(_run_mcal_one_chunk)(
            meds_files, start, end, seeds[chunk], mcal_config))

    with joblib.Parallel(
            n_jobs=n_chunks, backend='loky',
            verbose=50, max_nbytes=None) as p:
        outputs = p(jobs)

    assert not all([o is None for o in outputs]), (
        "All metacal fits failed!")

    output = eu.numpy_util.combine_arrlist(
        [o for o in outputs if o is not None])
    logger.info(' %d of %d metacal fits worked!', output.size, cat.size)

    mcal_pth = get_mcal_file_path(
        meds_dir=output_meds_dir,
        medsconf=MEDSCONF,
        tilename=tilename)
    logger.info(' metacal output: "%s"', mcal_pth)
    make_dirs_for_file(mcal_pth)
    fitsio.write(mcal_pth, output, clobber=True)

