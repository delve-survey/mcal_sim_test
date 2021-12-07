# copied directly from https://github.com/beckermr/misc/blob/68b42ea226913b97236ce5145378c6cf6fe05b8b/matts_misc/simple_des_y3_sims/band_infoing.py

import logging
import hashlib
import os

import numpy as np
import yaml
import desmeds

from constants import MEDSCONF, PIFF_RUN, POSITION_OFFSET, MAGZP_REF
# from .files import (
#    get_band_info_file,
#    make_dirs_for_file)
#from .des_info import add_extra_des_coadd_tile_info

logger = logging.getLogger(__name__)


def make_band_info(*, tilename, bands, output_meds_dir, n_files=None):
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
        'campaign': 'Y3A1_COADD',
        'source_type': 'finalcut',
        'piff_run': PIFF_RUN,
        'medsconf': MEDSCONF
    }
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


def make_dirs_for_file(filename):
    """Make all of the parent directories for a file at `filename`."""
    dirname = os.path.dirname(filename)
    if len(dirname) > 0:
        os.makedirs(dirname, exist_ok=True)


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


def get_piff_path_from_image_path(*, image_path, piff_run):
    """Get the piff path from the image path.
    Parameters
    ----------
    image_path : str
        A path to an immask file
    piff_run : str
        e.g. y3a1-v29
    Returns
    -------
    psf_path : str
        The path to the PSF.
    """
    img_bname = os.path.basename(image_path)
    piff_bname = img_bname.replace('immasked.fits.fz', 'piff.fits')
    expnum = int(piff_bname.split('_')[0][1:])

    exp_dir = os.path.join(
        '$PIFF_DATA_DIR',
        piff_run,
        str(expnum),
    )

    psf_path = os.path.join(
        exp_dir,
        piff_bname,
    )

    return psf_path


def add_extra_des_coadd_tile_info(*, info, piff_run):
    """Read the coadd tile info, load WCS info, and load PSF info for
    the DES Y3+ DESDM layout.
    Parameters
    ----------
    info: dict
        Info dict for a coadd tile
    piff_run : str
        The PIFF PSF run to use.
    Returns
    -------
    info : dict
        A dictionary with at least the following keys:
            'position_offset' : the offset to add to zero-indexed image
                coordinates to get transform them to the convention assumed
                by the WCS.
            'src_info' : list of dicts for the SE sources
            'image_path' : the path to the FITS file with the coadd image
            'image_ext' : the name of the FITS extension with the coadd image
            'weight_path' : the path to the FITS file with the coadd weight map
            'weight_ext' : the name of the FITS extension with the coadd weight
                map
            'bmask_path' : the path to the FITS file with the coadd bit mask
            'bmask_ext' : the name of the FITS extension with the coadd bit
                mask
            'seg_path' : the path to the FITS file with the coadd seg map
            'seg_ext' : the name of the FITS extension with the coadd seg map
            'image_flags' : any flags for the coadd image
            'scale' : a multiplicative factor to apply to the image
                (`*= scale`) and weight map (`/= scale**2`) for magnitude
                zero-point calibration.
            'magzp' : the magnitude zero point for the image
        The dictionaries in the 'src_info' list have at least the
        following keys:
            'image_path' : the path to the FITS file with the SE image
            'image_ext' : the name of the FITS extension with the SE image
            'bkg_path' : the path to the FITS file with the SE background image
            'bkg_ext' : the name of the FITS extension with the SE background
                image
            'weight_path' : the path to the FITS file with the SE weight map
            'weight_ext' : the name of the FITS extension with the SE weight
                map
            'bmask_path' : the path to the FITS file with the SE bit mask
            'bmask_ext' : the name of the FITS extension with the SE bit mask
            'psfex_path' : the path to the PSFEx PSF model
            'piff_path' : the path to the Piff PSF model
            'scale' : a multiplicative factor to apply to the image
                (`*= scale`) and weight map (`/= scale**2`) for magnitude
                zero-point calibration
            'magzp' : the magnitude zero point for the image
            'image_flags' : any flags for the SE image
            'position_offset' : the offset to add to zero-indexed image
                coordinates to get transform them to the convention assumed
                by the WCS.
    """

    info['position_offset'] = POSITION_OFFSET

    info['image_ext'] = 'sci'

    info['weight_path'] = info['image_path']
    info['weight_ext'] = 'wgt'

    info['bmask_path'] = info['image_path']
    info['bmask_ext'] = 'msk'

    info['seg_ext'] = 'sci'

    # always true for the coadd
    info['magzp'] = MAGZP_REF
    info['scale'] = 1.0
    info['image_shape'] = [10000, 10000]

    info['image_flags'] = 0

    for index, ii in enumerate(info['src_info']):
        ii['image_shape'] = [4096, 2048]
        ii['image_flags'] = 0

        ii['image_ext'] = 'sci'

        ii['weight_path'] = ii['image_path']
        ii['weight_ext'] = 'wgt'

        ii['bmask_path'] = ii['image_path']
        ii['bmask_ext'] = 'msk'

        ii['bkg_ext'] = 'sci'

        # wcs info
        ii['position_offset'] = POSITION_OFFSET

        # psfex psf
        ii['psfex_path'] = ii['psf_path']

        # piff
        ii['piff_path'] = get_piff_path_from_image_path(
            image_path=ii['image_path'],
            piff_run=piff_run,
        )

        # image scale
        ii['scale'] = 10.0**(0.4*(MAGZP_REF - ii['magzp']))


def get_piff_path_from_image_path(*, image_path, piff_run):
    """Get the piff path from the image path.
    Parameters
    ----------
    image_path : str
        A path to an immask file
    piff_run : str
        e.g. y3a1-v29
    Returns
    -------
    psf_path : str
        The path to the PSF.
    """
    img_bname = os.path.basename(image_path)
    piff_bname = img_bname.replace('immasked.fits.fz', 'piff.fits')
    expnum = int(piff_bname.split('_')[0][1:])

    exp_dir = os.path.join(
        '$PIFF_DATA_DIR',
        piff_run,
        str(expnum),
    )

    psf_path = os.path.join(
        exp_dir,
        piff_bname,
    )

    return psf_path
