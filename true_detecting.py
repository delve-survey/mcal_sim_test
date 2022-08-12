import os
import shutil
import logging

import yaml
import fitsio
import numpy as np

from files import (
    get_truth_catalog_path,
    get_band_info_file,
    expand_path,
    make_dirs_for_file)
from constants import MEDSCONF
from realistic_galaxying import init_descwl_catalog

logger = logging.getLogger(__name__)

TMP_DIR = os.environ['TMPDIR']

def make_true_detections(*, tilename, bands, output_meds_dir, box_size, config):
    """Make fake "true detection" catalogs as if source extractor had been
    run.

    Note this function also writes "fake" coadd image data so that
    MEDS files can be made.

    Parameters
    ----------
    tilename : str
        The DES coadd tile to run true detection.
    bands : str
        The bands to run true detection.
    output_meds_dir : str
        The output DESDATA/MEDS_DIR for the simulation data products.
    box_size : int
        The desired box size in the MEDS files made from the fake
        "true detection" catalogs. The source extractor columns are hacked
        so that the MEDS making code produces this box size.
    config : dict
        A dictionary with the config params for the simulations
    """

    logger.info(' processing coadd tile %s', tilename)

    tcat = fitsio.read(get_truth_catalog_path(
        meds_dir=output_meds_dir,
        medsconf=MEDSCONF,
        tilename=tilename))

    for band in bands:
        dest_cat_file = _copy_and_munge_coadd(
            tilename=tilename,
            band=band,
            output_meds_dir=output_meds_dir)

        _reformat_catalog(
            output_cat_file=dest_cat_file,
            truth_cat=tcat,
            tilename=tilename,
            band=band,
            box_size=box_size,
            config = config)


def _reformat_catalog(*, output_cat_file, truth_cat, tilename, band, box_size, config):
    
    # now we reformat to the sectractor output
    # note that we are hacking on the fields to force the
    # MEDS maker to use the right sized stamps
    # to do this we
    # 1. set the flux radius to zero
    # 2. set the x[y]min[max]_image fields to have the
    #   desired box size
    # 3. set b_world to 1 and a_world to 0
    # 4. set the flags to zero
        
    dtype = [
        ('number', 'i4'),
        ('xmin_image', 'i4'),
        ('ymin_image', 'i4'),
        ('xmax_image', 'i4'),
        ('ymax_image', 'i4'),
        ('x_image', 'f4'),
        ('y_image', 'f4'),
        ('alpha_j2000', 'f8'),
        ('delta_j2000', 'f8'),
        ('a_world', 'f4'),
        ('b_world', 'f4'),
        ('flags', 'i2'),
        ('flux_radius', 'f4')]
    srcext_cat = np.zeros(len(truth_cat), dtype=dtype)
    srcext_cat['number'] = np.arange(len(truth_cat)) + 1
    srcext_cat['x_image'] = truth_cat['x']
    srcext_cat['y_image'] = truth_cat['y']
    srcext_cat['alpha_j2000'] = truth_cat['ra']
    srcext_cat['delta_j2000'] = truth_cat['dec']
    srcext_cat['a_world'] = 1
    srcext_cat['b_world'] = 0
    srcext_cat['flags'] = 0
    srcext_cat['flux_radius'] = 0
    
    
    if config['gal_kws']['gal_source'] in ['simpleElliptical', 'varsize', 'varang', 'varsizeang', 'descwl']:
        
        srcext_cat['a_world'] = truth_cat['a_world']
        srcext_cat['b_world'] = truth_cat['b_world']
        srcext_cat['flux_radius'] = truth_cat['size']        


    half = int(box_size / 2)
    xint = (truth_cat['x'] + 0.5).astype(np.int32)
    srcext_cat['xmin_image'] = xint - half
    srcext_cat['xmax_image'] = (
        box_size - 1 + srcext_cat['xmin_image'])
    yint = (truth_cat['y'] + 0.5).astype(np.int32)
    srcext_cat['ymin_image'] = yint - half
    srcext_cat['ymax_image'] = (
        box_size - 1 + srcext_cat['ymin_image'])

    # now we add the new srcext catalog to the stash and
    # write it to disk
    make_dirs_for_file(output_cat_file)
    fitsio.write(output_cat_file, srcext_cat, clobber=True)


def _copy_and_munge_coadd(*, tilename, band, output_meds_dir):
    # read band info
    fname = get_band_info_file(
        meds_dir=output_meds_dir,
        medsconf=MEDSCONF,
        tilename=tilename,
        band=band)
    with open(fname, 'r') as fp:
        info = yaml.load(fp, Loader=yaml.Loader)


    dest_coadd_file = info['image_path'].replace(
        TMP_DIR, output_meds_dir)
    make_dirs_for_file(dest_coadd_file)
    dest_seg_file = info['seg_path'].replace(
        TMP_DIR, output_meds_dir)
    make_dirs_for_file(dest_seg_file)

    # copy over coadd file
    logger.info(' copying coadd and seg file for band %s', band)
    shutil.copy(
        expand_path(info['image_path']),
        dest_coadd_file)
    shutil.copy(
        expand_path(info['seg_path']),
        dest_seg_file)

    if dest_coadd_file.endswith('.fz'):
        logger.info(' decompressing coadd file for band %s', band)
        try:
            os.remove(dest_coadd_file[:-3])
        except Exception:
            pass

        os.system('funpack %s' % dest_coadd_file)

        try:
            os.remove(dest_coadd_file)
        except Exception:
            pass

        dest_coadd_file = dest_coadd_file[:-3]

    # write all zeros in the image
    logger.info(' zeroing coadd file for band %s', band)
    with fitsio.FITS(dest_coadd_file, mode='rw') as fp:
        fp[info['image_ext']].write(np.zeros((10000, 10000)))

    # repack
    logger.info(' compressing coadd file for band %s', band)
    os.system('fpack %s' % dest_coadd_file)

    try:
        os.remove(dest_coadd_file)
    except Exception:
        pass

    return info['cat_path'].replace(TMP_DIR, output_meds_dir)
