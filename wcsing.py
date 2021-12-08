from functools import lru_cache
import fitsio
import esutil.wcsutil
import galsim


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
