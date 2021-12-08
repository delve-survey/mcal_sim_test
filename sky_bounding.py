import numpy as np
from meds.bounds import Bounds
from meds.util import radec_to_uv


def get_rough_sky_bounds(
        *, im_shape, wcs, position_offset, bounds_buffer_uv, n_grid,
        celestial=True):
    """Get the rough boundry of a CCD on the sky for detecting if an object
    is on the CCD.

    Algorithm due to M. Jarvis w/ some changes from M. R. Becker.

    Parameters
    ----------
    im_shape : two-tuple of ints
        The shape of the image.
    wcs : `esutil.wcsutil.WCS` or `AffineWCS` object
        The wcs object that defines the transformation from pixels to the sky.
    position_offset : int
        The offset from zero-indexed pixels needed to use the WCS. For one-
        indexed pixels, pass 1.
    bounds_buffer_uv : float
        The buffer in arcseconds for the chip boundaries in (u, v) coordinates.
        A typical value is 16 arcseconds.
    n_grid : int
        Number of grid points to use in the small direction to construct
        the bounding box. A typical value is 4.
    celestial : bool, optional
        If True, interpret the WCS as a celestial WCS. Otherwise, treat it
        as an affine transformation. Default is True.

    Returns
    -------
    sky_bnds : `meds.bounds.Bounds`
        The bounding box in a spherical coordinate system centered on the CCD.
    ra_ccd : float
        The CCD center ra.
    dec_ccd : float
        The CCD center dec.

    Examples
    --------
    >>> sky_bnds, ra_ccd, dec_ccd = get_rough_sky_bounds(
    >>>     im_shape=(4096, 2048),
    >>>     wcs=wcs, position_offset=1, bounds_buffer_uv=16, n_grid=4)
    >>> # ra, dec are points to test
    >>> u, v = radec_to_uv(ra, dec, ra_ccd, dec_ccd)
    >>> in_sky_bnds = sky_bnds.contains_points(u, v)  # returs a bool mask
    >>> q = np.where(in_sky_bnds)
    """
    nrow, ncol = im_shape

    # set n_grid so that pixels are square-ish
    if ncol < nrow:
        n_grid_col = n_grid
        n_grid_row = np.ceil(float(nrow)/float(ncol))
    else:
        n_grid_row = n_grid
        n_grid_col = np.ceil(float(ncol)/float(nrow))

    # construct a grid
    # this is zero-indexed
    rows = np.arange(n_grid_row+1)*(nrow-1.0)/n_grid_row
    cols = np.arange(n_grid_col+1)*(ncol-1.0)/n_grid_col
    # we could save some time by just doing the outside but shrug
    rows, cols = np.meshgrid(rows, cols)
    rows = rows.ravel()
    cols = cols.ravel()

    # get ra, dec
    ra, dec = wcs.image2sky(
        x=cols + position_offset,
        y=rows + position_offset)

    # get ccd center
    # these are zero-indexed
    row_ccd = (nrow - 1)/2
    col_ccd = (ncol - 1)/2
    ra_ccd, dec_ccd = wcs.image2sky(
        x=col_ccd + position_offset,
        y=row_ccd + position_offset)

    if celestial:
        # get u,v - ccd is at 0,0 by def
        u, v = radec_to_uv(ra, dec, ra_ccd, dec_ccd)

        # build bounds with buffer and cos(dec) factors
        vrad = np.deg2rad(v / 3600.0)  # arcsec to degrees
        ufac = np.cos(vrad).min()

        ubuff = bounds_buffer_uv / ufac
        vbuff = bounds_buffer_uv
    else:
        u = ra - ra_ccd
        v = dec - dec_ccd
        ubuff = bounds_buffer_uv
        vbuff = bounds_buffer_uv

    sky_bnds = Bounds(u.min() - ubuff,
                      u.max() + ubuff,
                      v.min() - vbuff,
                      v.max() + vbuff)

    return sky_bnds, ra_ccd, dec_ccd
