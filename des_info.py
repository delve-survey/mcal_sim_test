from constants import MAGZP_REF, POSITION_OFFSET
from files import get_piff_path_from_image_path


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
        if piff_run is not None:
            ii['piff_path'] = get_piff_path_from_image_path(
                image_path=ii['image_path'],
                piff_run=piff_run,
            )

        # image scale
        ii['scale'] = 10.0**(0.4*(MAGZP_REF - ii['magzp']))
