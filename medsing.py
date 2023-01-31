import json
import os
import logging
import tempfile

import numpy as np
import meds.util
import fitsio
import galsim
import yaml
from esutil.ostools import StagedOutFile
from meds.maker import MEDSMaker
import desmeds.util

from constants import MAGZP_REF, MEDSCONF
from psf_wrapper import PSFWrapper
from wcsing import get_galsim_wcs
from files import (
    get_band_info_file, get_meds_file_path, expand_path, make_dirs_for_file)

logger = logging.getLogger(__name__)


TMP_DIR = os.environ['TMPDIR']

def make_meds_files(*, tilename, bands, output_meds_dir, psf_kws, meds_config):
    """Make a MEDS file for a given band and tilename.

    Parameters
    ----------
    tilename : str
        The DES coadd tile to run true detection.
    bands : str
        The bands to run true detection.
    output_meds_dir : str
        The output DEADATA/MEDS_DIR for the simulation data products.
    psf_kws : dict
        The dictionary of PSF config information.
    meds_config : dict
        The MEDS making configuration file. See the default one in
        `work/simple_des_y3_sims/default_configs/meds.yaml`.
    """

    logger.info(' making meds files for coadd tile %s', tilename)

    # force this
    meds_config['magzp_ref'] = MAGZP_REF
    meds_config['psf'] = {'type': 'psfex'}
    meds_config['use_joblib'] = True

    # read info files
    info = {}
    for band in bands:
        # get info about files
        fname = get_band_info_file(
            meds_dir=output_meds_dir,
            medsconf=MEDSCONF,
            tilename=tilename,
            band=band)
        with open(fname, 'r') as fp:
            info[band] = yaml.load(fp, Loader=yaml.Loader)
            print(info[band]['image_path'])
            #for k in info[band].keys():
            #    print(k, info[band][k])
    # always get the truth catalog from r band
    cat = fitsio.read(info['r']['cat_path'].replace(
        TMP_DIR, output_meds_dir))

    for band in bands:
        logger.info(' doing band %s', band)

        # get all of the components for the file
        obj_data = _make_meds_input_data_struct(
            cat=cat,
            allowed_box_sizes=meds_config['allowed_box_sizes'],
            min_box_size=meds_config['min_box_size'],
            max_box_size=meds_config['max_box_size'],
            sigma_fac=meds_config['sigma_fac'])
        image_info = _make_meds_image_info_struct(
            info=info[band], output_meds_dir=output_meds_dir)
        meta_data = _make_meds_metadata(band=band, tilename=tilename)
        psf_data = _build_psf_data(
            info=info[band],
            psf_kws=psf_kws,
            output_meds_dir=output_meds_dir)

        # make the file in a tmp dir and then stage out
        maker = MEDSMaker(
            obj_data,
            image_info,
            psf_data=psf_data,
            config=meds_config,
            meta_data=meta_data)

        final_meds_file = get_meds_file_path(
            meds_dir=output_meds_dir,
            medsconf=MEDSCONF,
            tilename=tilename,
            band=band)
        make_dirs_for_file(final_meds_file)

        with tempfile.TemporaryDirectory() as tmpdir:
            with StagedOutFile(final_meds_file, tmpdir=tmpdir) as sf:
                uncompressed_file = sf.path.replace('.fits.fz', '.fits')
                make_dirs_for_file(uncompressed_file)
                maker.write(uncompressed_file)

                # make sure to remove the destination file when fpacking
                try:
                    os.remove(sf.path)
                except Exception:
                    pass
                desmeds.util.fpack_file(uncompressed_file)
                try:
                    os.remove(uncompressed_file)
                except Exception:
                    pass


def _build_psf_data(*, info, psf_kws, output_meds_dir):
    def _load_psf_data(_info, force_gauss=False):
        wcs = get_galsim_wcs(
            image_path=_info['image_path'].replace(
                TMP_DIR, output_meds_dir),
            image_ext=_info['image_ext'])
        
        if psf_kws['type'] == 'gauss' or force_gauss:
            return PSFWrapper(galsim.Gaussian(fwhm=0.9), wcs)
        
        #elif psf_kws['type'] == 'piff':
        #    from ..des_piff import DES_Piff
        #    piff_model = DES_Piff(expand_path(_info['piff_path']))
        #    return PSFWrapper(piff_model, wcs)
        
        elif psf_kws['type'] == 'gauss-pix':
            from gauss_pix_psf import GaussPixPSF
            kwargs = {k: psf_kws[k] for k in psf_kws if k != 'type'}
            psf_model = GaussPixPSF(**kwargs)
            return PSFWrapper(psf_model, wcs)

        elif psf_kws['type'] == 'nongauss-pix':
            from nongauss_pix_psf import NonGaussPixPSF
            kwargs = {k: psf_kws[k] for k in psf_kws if k != 'type'}
            psf_model = NonGaussPixPSF(**kwargs)
            return PSFWrapper(psf_model, wcs)
        
        elif psf_kws['type'] == 'psfex':
            from galsim.des import DES_PSFEx
            psfex_model = DES_PSFEx(expand_path(_info['psfex_path']), wcs = wcs)
            return PSFWrapper(psfex_model, wcs)
        
        elif psf_kws['type'] == 'des_psfex':
            from des_psfex import DES_PSFEx_Deconv
            psfex_model = DES_PSFEx_Deconv(expand_path(_info['psfex_path']), wcs = wcs)
            return PSFWrapper(psfex_model, wcs)
        
        
        elif psf_kws['type'] == 'psfex_deconvolved':
            from psfex_deconvolved import PSFEx_Deconv
            psfex_model = PSFEx_Deconv(expand_path(_info['psfex_path']), wcs = wcs)
            return PSFWrapper(psfex_model, wcs)
        
        else:
            raise ValueError("psf type '%s' is not valid!" % psf_kws['type'])

    force_gauss = psf_kws['type'] in ['psfex', 'psfex_deconvolved', 'des_psfex', 'piff']
    psf_data = [_load_psf_data(info, force_gauss=force_gauss)]
    for se_info in info['src_info']:
        #print(se_info.keys())
        psf_data.append(_load_psf_data(se_info))
    return psf_data


def _make_meds_metadata(*, band, tilename):
    meta = np.zeros(1, dtype=[
        ('magzp_ref', 'f8'),
        ('band', 'S1'),
        ('tilename', 'S12')])
    meta['magzp_ref'] = MAGZP_REF
    meta['band'] = band
    meta['tilename'] = tilename
    return meta


def _make_meds_image_info_struct(*, info, output_meds_dir):
    def _munge_path(pth):
        return pth.replace(TMP_DIR, output_meds_dir)

    # get WCS structures
    wcs_json = _load_wcs_json(info=info, output_meds_dir=output_meds_dir)
    wcs_len = max([len(j) for j in wcs_json])

    # compute the max path length
    path_len = [
        len(_munge_path(info['image_path'])),
        len(_munge_path(info['weight_path'])),
        len(_munge_path(info['seg_path'])),
        len(_munge_path(info['bmask_path']))]
    for se_info in info['src_info']:
        path_len += [
            len(_munge_path(se_info['image_path'])),
            len(_munge_path(se_info['weight_path'])),
            len(_munge_path(se_info['bkg_path'])),
            len(_munge_path(se_info['bmask_path']))]
    path_len = max(path_len)

    # now fill the array
    dtype = meds.util.get_image_info_dtype(
        path_len,
        wcs_len=wcs_len,
        ext_len=3)
    image_info = np.zeros(len(info['src_info']) + 1, dtype=dtype)

    image_info['image_id'] = np.arange(len(image_info))
    image_info['image_flags'] = 0
    image_info['position_offset'] = 1
    for i, wj in enumerate(wcs_json):
        image_info['wcs'][i] = wj

    image_info['scale'][0] = info['scale']
    image_info['magzp'][0] = info['magzp']
    image_info['image_path'][0] = _munge_path(info['image_path'])
    image_info['image_ext'][0] = info['image_ext']
    image_info['weight_path'][0] = _munge_path(info['weight_path'])
    image_info['weight_ext'][0] = info['weight_ext']
    image_info['bmask_path'][0] = _munge_path(info['bmask_path'])
    image_info['bmask_ext'][0] = info['bmask_ext']
    image_info['bkg_path'][0] = ""
    image_info['bkg_ext'][0] = ""
    image_info['seg_path'][0] = _munge_path(info['seg_path'])
    image_info['seg_ext'][0] = info['seg_ext']
    for i, se_info in enumerate(info['src_info']):
        image_info['scale'][i+1] = se_info['scale']
        image_info['magzp'][i+1] = se_info['magzp']
        image_info['image_path'][i+1] = _munge_path(se_info['image_path'])
        image_info['image_ext'][i+1] = se_info['image_ext']
        image_info['weight_path'][i+1] = _munge_path(se_info['weight_path'])
        image_info['weight_ext'][i+1] = se_info['weight_ext']
        image_info['bmask_path'][i+1] = _munge_path(se_info['bmask_path'])
        image_info['bmask_ext'][i+1] = se_info['bmask_ext']
        image_info['bkg_path'][i+1] = _munge_path(se_info['bkg_path'])
        image_info['bkg_ext'][i+1] = se_info['bkg_ext']
        image_info['seg_path'][i+1] = ""
        image_info['seg_ext'][i+1] = ""

    return image_info


def _load_wcs_json(*, info, output_meds_dir):
    def _munge_header(hd):
        return {k.lower(): hd[k] for k in hd if k is not None}

    wcs_json = []
    hd = fitsio.read_header(
        info['image_path'].replace(TMP_DIR, output_meds_dir),
        ext=info['image_ext'])
    wcs_json.append(json.dumps(_munge_header(hd)))
    for se_info in info['src_info']:
        hd = fitsio.read_header(
            se_info['image_path'].replace(TMP_DIR, output_meds_dir),
            ext=se_info['image_ext'])
        wcs_json.append(json.dumps(_munge_header(hd)))

    return wcs_json


def _make_meds_input_data_struct(
        *, cat, allowed_box_sizes, min_box_size, max_box_size, sigma_fac):
    """Make the input data structure for the MEDS maker.

    Parameters
    ----------
    cat : np.ndarray
        The coadd catalog.
    allowed_box_sizes : list of ints
        A list of the allowed postage stamp box sizes.
    min_box_size : int
        The minimum allowed box size. This value overrides any smaller values
        in `allowed_box_sizes`.
    max_box_size : int
        The maximum allowed box size. This value overrides any larger values
        in `allowed_box_sizes`.
    sigma_fac : float
        The factor by which to scale the flux radius. A value around 5 is
        standard.

    Returns
    -------
    input_data : np.ndarray
        The input structured array for the MEDS maker.
    """
    dtype = meds.util.get_meds_input_dtype(extra_fields=[('number', 'i8')])
    input_data = np.zeros(len(cat), dtype=dtype)

    try:
        input_data['id'] = cat['number']
        input_data['number'] = cat['number']
        input_data['ra'] = cat['alpha_j2000']
        input_data['dec'] = cat['delta_j2000']
    except:
        input_data['id'] = cat['number'.upper()]
        input_data['number'] = cat['number'.upper()]
        input_data['ra'] = cat['alpha_j2000'.upper()]
        input_data['dec'] = cat['delta_j2000'.upper()]
 
    input_data['box_size'] = _get_box_sizes(
        cat=cat,
        allowed_box_sizes=allowed_box_sizes,
        min_box_size=min_box_size,
        max_box_size=max_box_size,
        sigma_fac=sigma_fac)
    return input_data


def _get_box_sizes(
        *, cat, allowed_box_sizes, min_box_size, max_box_size, sigma_fac):
    """Get the box sizes for the coadd catalog.

    Parameters
    ----------
    cat : np.ndarray
        The coadd catalog.
    allowed_box_sizes : list of ints
        A list of the allowed postage stamp box sizes.
    min_box_size : int
        The minimum allowed box size. This value overrides any smaller values
        in `allowed_box_sizes`.
    max_box_size : int
        The maximum allowed box size. This value overrides any larger values
        in `allowed_box_sizes`.
    sigma_fac : float
        The factor by which to scale the flux radius. A value around 5 is
        standard.

    Returns
    -------
    box_sizes : np.ndarray
        The array of box sizes.
    """
    sigma_size = get_sigma_size(cat=cat, sigma_fac=sigma_fac)

    # now do row and col sizes
    try:
        row_size = cat['ymax_image'] - cat['ymin_image'] + 1
        col_size = cat['xmax_image'] - cat['xmin_image'] + 1
    except:
        row_size = cat['ymax_image'.upper()] - cat['ymin_image'.upper()] + 1
        col_size = cat['xmax_image'.upper()] - cat['xmin_image'.upper()] + 1

    # get max of all three
    box_size = np.vstack(
        (col_size, row_size, sigma_size)).max(axis=0)

    # clip to range
    box_size = box_size.clip(min_box_size, max_box_size)

    # now put in fft sizes
    bins = [0]
    bins.extend([sze for sze in allowed_box_sizes
                 if sze >= min_box_size and sze <= max_box_size])

    if bins[-1] != max_box_size:
        bins.append(max_box_size)

    bin_inds = np.digitize(box_size, bins, right=True)
    bins = np.array(bins)

    return bins[bin_inds]


def get_sigma_size(*, cat, sigma_fac):
    """Get an object size based on its flux radius and ellipticity.

    Parameters
    ----------
    cat : np.ndarray
        The coadd catalog.
    sigma_fac : float
        The factor by which to scale the flux radius. A value around 5 is
        standard.

    Returns
    -------
    sigma_size : np.ndarray
        The array of sizes.
    """

    fwhm_fac = 2*np.sqrt(2*np.log(2))
    
    try:
        ellipticity = 1.0 - cat['b_world'] / cat['a_world']
        sigma = cat['flux_radius'] * 2.0 / fwhm_fac
    except:
        ellipticity = 1.0 - cat['b_world'.upper()] / cat['a_world'.upper()]
        sigma = cat['flux_radius'.upper()] * 2.0 / fwhm_fac 
    drad = sigma * sigma_fac
    drad = drad * (1.0 + ellipticity)
    drad = np.ceil(drad)
    # sigma size is twice the radius
    sigma_size = 2 * drad.astype('i4')

    return sigma_size
