import os
import functools
import collections
import numpy as np

import galsim
import fitsio

# @functools.lru_cache(maxsize=8)
def _cached_catalog_read():
    stars    = os.path.join(os.environ.get('LSSTSTARSIM_DIR', '.'), 'LSST_Stars.fits',)
    binaries = os.path.join(os.environ.get('LSSTSTARSIM_DIR', '.'), 'LSST_Binaries.fits',)
    
    return fitsio.read(stars), fitsio.read(binaries)


def init_lsst_starsim_catalog(*, rng):
    
    stars, binaries = _cached_catalog_read()
    
    return stars, binaries