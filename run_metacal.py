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


def _run_mcal_one_chunk(meds_files, start, end, seed, mcal_config):
    """Run metcal for `meds_files` only for objects from `start` to `end`.

    Note that `start` and `end` follow normal python indexing conventions so
    that the list of indices processed is `list(range(start, end))`.

    Parameters
    ----------
    meds_files : list of str
        A list of paths to the MEDS files.
    start : int
        The starting index of objects in the file on which to run metacal.
    end : int
        One plus the last index to process.
    seed : int
        The seed for the RNG.
    mcal_config : yaml
        The config file for the metacal run

    Returns
    -------
    output : np.ndarray
        The metacal outputs.
    """
    rng = np.random.RandomState(seed=seed)

    # seed the global RNG to try to make things reproducible
    np.random.seed(seed=rng.randint(low=1, high=2**30))

    output = None
    mfiles = []
    data = []
    try:
        # get the MEDS interface
        for m in meds_files:
            mfiles.append(NGMixMEDS(m))
        mbmeds = MultiBandNGMixMEDS(mfiles)
        cat = mfiles[0].get_cat()

        for ind in range(start, end):
            o = mbmeds.get_mbobs(ind)
            
            o = _strip_coadd(o, mcal_config) #Remove coadd since it isnt used in fitting
            o = _strip_zero_flux(o, mcal_config) #Remove any obs with zero flux
            
            if mcal_config['custom']['maxbadfrac'] > 0: o = _strip_10percent_masked(o, mcal_config) #Remove obs with many bad pixs
            if mcal_config['custom']['goodfrac']: o = _get_masked_frac(o, mcal_config) #Get gauss-weighted fraction of good pixels
            if mcal_config['custom']['symmetrize_mask']: o = _symmetrize_mask(o, mcal_config) #Symmetrize bmask for 180 deg symmetry
            if mcal_config['custom']['interp_bad_pixels']: o = _fill_empty_pix(o, rng, mcal_config) #Interpolate empty pixels
            
            o = _apply_pixel_scale(o, mcal_config) #Not sure??

            skip_me = False
            for ol in o:
                if len(ol) == 0:
                    logger.debug(' not all bands have images - skipping!')
                    skip_me = True
            if skip_me:
                continue

            o.meta['id'] = ind
            o[0].meta['Tsky'] = 1
            o[0].meta['magzp_ref']   = MAGZP_REF
            o[0][0].meta['orig_col'] = cat['orig_col'][ind, 0]
            o[0][0].meta['orig_row'] = cat['orig_row'][ind, 0]

            if mcal_config['custom']['goodfrac']:
                #put all the good_fraction numbers into one list
                #one entry per cutout (so all bands are combined here)
                good_frac = []
                weight    = []
                for _one in o:
                    for _two in _one:
                        good_frac.append(_two.meta['good_frac'])
                        weight.append(_two.meta['weight'])
                    
            nband = len(o)
            mcal = MetacalFitter(mcal_config, nband, rng)

            try:
                mcal.go([o])
                res = mcal.result
            except GMixRangeError as e:
                logger.debug(" metacal error: %s", str(e))
                res = None

            if res is not None:
                if mcal_config['custom']['goodfrac']:
                    res['good_frac'] = np.average(good_frac, weights = weight) #store mean good_fraction per object
                else:
                    res['good_frac'] = 1 #Else, we're assuming image is "perfect" (completely unmasked) == 1
                data.append(res)

        if len(data) > 0:
            output = eu.numpy_util.combine_arrlist(data)
    finally:
        for m in mfiles:
            m.close()

    return output


def _strip_coadd(mbobs, mcal_config):
    _mbobs = MultiBandObsList()
    _mbobs.update_meta_data(mbobs.meta)
    for ol in mbobs:
        _ol = ObsList()
        _ol.update_meta_data(ol.meta)
        for i in range(1, len(ol)):
            _ol.append(ol[i])
        _mbobs.append(_ol)
    return _mbobs


def _strip_zero_flux(mbobs, mcal_config):
    _mbobs = MultiBandObsList()
    _mbobs.update_meta_data(mbobs.meta)
    for ol in mbobs:
        _ol = ObsList()
        _ol.update_meta_data(ol.meta)
        for i in range(len(ol)):
            if np.sum(ol[i].image) > 0:
                _ol.append(ol[i])
        _mbobs.append(_ol)
    return _mbobs


def _apply_pixel_scale(mbobs, mcal_config):
    for ol in mbobs:
        for o in ol:
            scale    = o.jacobian.get_scale()
            scale2   = scale * scale
            scale4   = scale2 * scale2
            o.image  = o.image / scale2
            o.weight = o.weight * scale4
            
            if mcal_config['custom']['interp_bad_pixels']:
                o.noise  = o.noise / scale2
            
    return mbobs

def _strip_10percent_masked(mbobs, mcal_config):
    _mbobs = MultiBandObsList()
    _mbobs.update_meta_data(mbobs.meta)
    
    #Loop over different band observations (r, i, z)
    for ol in mbobs:
        _ol = ObsList()
        _ol.update_meta_data(ol.meta)
        
        #Loop over different exposures/cutouts in each band
        for i in range(len(ol)):
            
            msk = ol[i].bmask.astype(bool) #Mask where TRUE means bad pixel
            
#             print("strip", np.average(msk), mcal_config['custom']['maxbadfrac'])
                
            if np.average(msk) >= mcal_config['custom']['maxbadfrac']:
                continue
            
            _ol.append(ol[i])
        _mbobs.append(_ol)
    return _mbobs

def _get_masked_frac(mbobs, mcal_config):
    
    _mbobs = MultiBandObsList()
    _mbobs.update_meta_data(mbobs.meta)
    
    gauss = galsim.Gaussian(fwhm = 1.2) #Fixed aperture gauss weights for image
    
    #Loop over different band observations (r, i, z)
    for ol in mbobs:
        _ol = ObsList()
        _ol.update_meta_data(ol.meta)
        
        #Loop over different exposures/cutouts in each band
        for i in range(len(ol)):
            
            msk = ol[i].bmask.astype(bool) #Mask where TRUE means bad pixel
            wgt = np.median(ol[i].weight[ol[i].weight != 0]) #Median weight used to populate noise in empty pix
            
            #get wcs of this observations
            wcs = ol[i].jacobian.get_galsim_wcs()

            #Create gaussian weights image (as array)
            gauss_wgt = gauss.drawImage(nx = msk.shape[0], ny = msk.shape[1], wcs = wcs, method = 'real_space').array 

            #msk is nonzero for bad pixs. Invert it, and convert to int
            good_frac = np.average(np.invert(msk).astype(int), weights = gauss_wgt) #Fraction of missing values

            #Save fraction of good pix. Will use later to remove
            #problematic objects directly from metacal catalog
            ol[i].meta['good_frac'] = good_frac
            ol[i].meta['weight']    = wgt
            
#             print("goodfrac", good_frac, wgt)

            _ol.append(ol[i])
        _mbobs.append(_ol)
    return _mbobs

def _symmetrize_mask(mbobs, mcal_config):
    
    _mbobs = MultiBandObsList()
    _mbobs.update_meta_data(mbobs.meta)
    
    #Loop over different band observations (r, i, z)
    for ol in mbobs:
        _ol = ObsList()
        _ol.update_meta_data(ol.meta)
        
        #Loop over different exposures/cutouts in each band
        for i in range(len(ol)):
            
            msk = ol[i].bmask.astype(bool) #Mask where TRUE means bad pixel
                
#             print("symm1", msk.sum(), np.sum(ol[i].bmask))
            #Rotate because Metacal needs this
            msk |= np.rot90(msk, k = 1)
            
            #Write rotated mask back to observation
            ol[i].bmask = msk.astype(np.int32)
            
#             print("symm2", msk.sum(), np.sum(ol[i].bmask))
            
            _ol.append(ol[i])
        _mbobs.append(_ol)
    return _mbobs
    
def _fill_empty_pix(mbobs, rng, mcal_config):
    _mbobs = MultiBandObsList()
    _mbobs.update_meta_data(mbobs.meta)
    
    #Loop over different band observations (r, i, z)
    for ol in mbobs:
        _ol = ObsList()
        _ol.update_meta_data(ol.meta)
        
        #Loop over different exposures/cutouts in each band
        for i in range(len(ol)):
            
            msk = ol[i].bmask.astype(bool) #Mask where TRUE means bad pixel
            wgt = np.median(ol[i].weight[ol[i].weight != 0]) #Median weight used to populate noise in empty pix
            
            #Observation doesn't have noise image, and so add noise image in.
            #Just random gaussian noise image using weights
            #Need to do this for interpolation step    
            ol[i].noise = rng.normal(loc = 0, scale = 1/np.sqrt(wgt), size = ol[i].image.shape)

            
            #If there are any bad mask pixels, then do interpolation
            if np.any(msk):
                
                #Interpolate image to fill in gaps. Setting maxfrac=1 since maxfrac is checked beforehand
                im    = interpolate_image_at_mask(image=ol[i].image, weight=wgt, bad_msk=msk, 
                                                  rng=rng, maxfrac=1, buff=4,
                                                  fill_isolated_with_noise=True)

                #Interpolate over noise image
                noise = interpolate_image_at_mask(image=ol[i].noise, weight=wgt, bad_msk=msk, 
                                                  rng=rng, maxfrac=1, buff=4,
                                                  fill_isolated_with_noise=True)

                #If we can't interpolate image or noise due to lack of data
                #then we skip this observation (it is stripped from MultiBandObs list)
                if (im is None) | (noise is None):
                    continue
                    
                #Set all masked pixel weights to 0.0
                ol[i].image  = im
                ol[i].weight = np.where(msk, 0, ol[i].weight)
                ol[i].noise  = noise

            
            _ol.append(ol[i])
        _mbobs.append(_ol)
    return _mbobs