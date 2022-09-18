import os
from contextlib import contextmanager


@contextmanager
def environ(**env):
    """Temporarily set environment variables inside the context manager and
    fully restore previous environment afterwardsself.

    Pass the env vars as keywords.

    Code from https://gist.github.com/igniteflow/7267431#gistcomment-2553451.

    Example
    -------
    >>> with environ(USER='xyz'):
    >>>     print(os.environ['USER'])
    >>> print(os.environ['USER'])
    """
    original_env = {key: os.getenv(key) for key in env}
    os.environ.update(env)
    try:
        yield
    finally:
        for key, value in original_env.items():
            if value is None:
                del os.environ[key]
            else:
                os.environ[key] = value


def expand_path(path):
    """Perform a full unix-like path expansion on `path`."""
    return os.path.expandvars(os.path.expanduser(path))


def make_dirs_for_file(filename):
    """Make all of the parent directories for a file at `filename`."""
    dirname = os.path.dirname(filename)
    if len(dirname) > 0:
        os.makedirs(dirname, exist_ok=True)


def get_mcal_file_path(*, meds_dir, medsconf, tilename):
    """Get the path to the output metacal file.

    Parameters
    ----------
    meds_dir : str
        The DESDATA/MEDS_DIR path where the info file is located.
    medsconf : str
        The MEDS file version (e.g., 'y3v02').
    tilename : str
        The DES coadd tilename (e.g., 'DES2122+0001').

    Returns
    -------
    mcal_pth : str
        The path to the metacal file.
    """
    return os.path.join(
        meds_dir,
        'metacal',
        medsconf,
        "%s_metacal.fits" % tilename)


def get_meds_file_path(*, meds_dir, medsconf, tilename, band):
    """Get the MEDS file for a given tile and band.

    Parameters
    ----------
    meds_dir : str
        The DESDATA/MEDS_DIR path where the info file is located.
    medsconf : str
        The MEDS file version (e.g., 'y3v02').
    tilename : str
        The DES coadd tilename (e.g., 'DES2122+0001').
    band : str
        The band as a string (e.g., `'r'`).

    Returns
    -------
    meds_file : str
        The path to the MEDS file.
    """
    return os.path.join(
        meds_dir,
        'meds',
        medsconf,
        tilename,
        "%s_%s_meds-%s.fits.fz" % (tilename, band, medsconf))


def get_truth_catalog_path(*, meds_dir, medsconf, tilename):
    """Get the truth catalog path.

    Parameters
    ----------
    meds_dir : str
        The DESDATA/MEDS_DIR path where the info file is located.
    medsconf : str
        The MEDS file version (e.g., 'y3v02').
    tilename : str
        The DES coadd tilename (e.g., 'DES2122+0001').

    Returns
    -------
    truth_file_path : str
        The path to the truth file.
    """
    return os.path.join(
        meds_dir,
        'simple_des_y3_sims',
        medsconf,
        'truthcats',
        '%s_truthcat.fits' % tilename
    )


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


def get_swarp_files_path(*, meds_dir, medsconf):
    """Get the path of files related to swarp (txt files of exp lists
    and then also the final detection coadd),
    that go into the swarp coadd for this `tilename` and `band`.

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
    path : str
        The path of list files for coadd generation.
    """
    return os.path.join(
        meds_dir,
        'simple_des_y3_sims',
        medsconf,
        'swarp_files')


def get_nwgint_path(*, meds_dir, medsconf, band):
    """Get the path of the ngwint files used as intermediaries
    that go into the swarp coadd for this `tilename` and `band`.

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
    path : str
        The path of list files for coadd generation.
    """
    return os.path.join(
        meds_dir,
        'simple_des_y3_sims',
        medsconf,
        'nwgint',
        'sources-%s'%band)


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
