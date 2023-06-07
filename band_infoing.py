import logging
import hashlib

import numpy as np
import yaml
import desmeds

from constants import MEDSCONF, PIFF_RUN
from files import (
    get_band_info_file,
    make_dirs_for_file)
from des_info import add_extra_des_coadd_tile_info

logger = logging.getLogger(__name__)


def make_band_info(*, tilename, bands, output_meds_dir, config, n_files=None):
    """Make YAML files with the information on each band.

    Parameters
    ----------
    tilename : str
        The DES coadd tilename (e.g., 'DES2122+0001').
    bands : list of str
        A list of bands to process (e.g., `['r', 'i', 'z']`).
    output_meds_dir : str
        The DESDATA/MEDS_DIR path where the info file should be written.
    n_files : int, optional
        If not `None`, then only keep this many files for the sources. Useful
        for testing.

    Returns
    -------
    info : dict
        A dictionary mapping the band name to the info file. Note that these
        paths are relative to the environment variable '$MEDS_DIR' in the
        returned file path. Replace this with `output_meds_dir` to read
        the file.
    """

    logger.info(' processing coadd tile %s', tilename)

    cfg = {
        'source_type': 'finalcut',
        'piff_run': PIFF_RUN,
        'medsconf': MEDSCONF
        }
    
    if config['name'] == 'DES':
        cfg['campaign'] = 'Y3A1_COADD'
    elif config['name'] == 'DECADE':
        cfg['campaign'] = 'DR3_1'
    
    fnames = {}
    for band in bands:
        band_info_file = get_band_info_file(
            meds_dir=output_meds_dir,
            medsconf=cfg['medsconf'],
            tilename=tilename,
            band=band)
        prep = desmeds.desdm_maker.Preparator(
                cfg,
                tilename,
                band,
            )
        prep.coadd["piff_campaign"] = None
        if prep.coadd.sources is not None:
            prep.coadd.sources["piff_campaign"] = None
        prep.go()
        info = prep.coadd.get_info()
        add_extra_des_coadd_tile_info(info=info, piff_run=cfg['piff_run'])

        # build hashes and sort
        hashes = []
        for i in range(len(info['src_info'])):
            hash_str = "%s%s" % (
                info['src_info'][i]['expnum'],
                info['src_info'][i]['ccdnum'])
            hashes.append(hashlib.md5(hash_str.encode('utf-8')).hexdigest())
        inds = np.argsort(hashes)
        new_src_info = [info['src_info'][i] for i in inds]
        info['src_info'] = new_src_info

        if n_files is not None:
            info['src_info'] = info['src_info'][:n_files]

        make_dirs_for_file(band_info_file)
        with open(band_info_file, 'w') as fp:
            yaml.dump(info, fp)

        fnames[band] = band_info_file.replace(output_meds_dir, '$MEDS_DIR')

    return fnames
