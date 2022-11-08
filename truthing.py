import numpy as np


def make_coadd_grid_radec(*, n_grid, coadd_wcs, rng, return_xy=False):
    """Make a grid of points in the coadd image coordinate system and
    return their locations in ra-dec.

    Parameters
    ----------
    n_grid : int
        The number of objects across the grid in each direction. The total
        number of objects will be `n_grid**2`.
    coadd_wcs : esutil.wcsutil.WCS
        The coadd WCS solution.
    rng : np.random.RandomState
        An RNG to use. This RNg is used to dither the locations on the coadd
        grid within a pixel.
    return_xy : bool, optional
        If True, also return the x and y positions. Default is False

    Returns
    -------
    ra : np.ndarray
        The array of ra positions of the sources.
    dec : np.ndarray
        The array of dec positions of the sources.
    x : np.ndarray
        The array of column positions. Only returned if `return_xy=True`.
    y : np.ndarray
        The array of row positions. Only returned if `return_xy=True`.
    """
    L = 10000  # hard code this since it will not change
    dL = L / n_grid
    dL_2 = dL / 2

    x = []
    y = []
    for row_ind in range(n_grid):
        for col_ind in range(n_grid):
            _x = col_ind * dL + dL_2 + 1
            _y = row_ind * dL + dL_2 + 1

            # dither
            _x += rng.uniform(low=-0.5, high=0.5)
            _y += rng.uniform(low=-0.5, high=0.5)

            x.append(_x)
            y.append(_y)

    x = np.array(x)
    y = np.array(y)
    ra, dec = coadd_wcs.image2sky(x, y)

    if return_xy:
        return ra, dec, x, y
    else:
        return ra, dec
    
    
def make_coadd_random_radec(*, n_gal, coadd_wcs, rng, return_xy=False):
    """Make a set of random points in the coadd image coordinate system and
    return their locations in ra-dec.

    Parameters
    ----------
    n_gal : int
        The number of objects across the grid. The total
        number of objects will be `n_gal`.
    coadd_wcs : esutil.wcsutil.WCS
        The coadd WCS solution.
    rng : np.random.RandomState
        An RNG to use. This RNg is used to dither the locations on the coadd
        grid within a pixel.
    return_xy : bool, optional
        If True, also return the x and y positions. Default is False

    Returns
    -------
    ra : np.ndarray
        The array of ra positions of the sources.
    dec : np.ndarray
        The array of dec positions of the sources.
    x : np.ndarray
        The array of column positions. Only returned if `return_xy=True`.
    y : np.ndarray
        The array of row positions. Only returned if `return_xy=True`.
    """
    L = 10000  # hard code this since it will not change

    x = rng.uniform(low = 0, high = L, size = n_gal)
    y = rng.uniform(low = 0, high = L, size = n_gal)
    
    #No dither because positions are already random 
    #x += rng.uniform(low=-0.5, high=0.5)
    #y += rng.uniform(low=-0.5, high=0.5)

    ra, dec = coadd_wcs.image2sky(x, y)

    if return_xy:
        return ra, dec, x, y
    else:
        return ra, dec
