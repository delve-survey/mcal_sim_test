import numpy as np
import galsim
import joblib


def render_sources_for_image(
        *, image_shape, wcs, draw_method, src_inds, src_func, n_jobs=None):
    """Render a list of sources for a single image.

    Parameters
    ----------
    image_shape : tuple of ints
        The shape of the final image.
    wcs : galsim WCS object
        The WCS for the image.
    draw_method : str
        The method used to draw the image. See the docs of `GSObject.drawImage`
        for details and options. Usually 'auto' is correct unless using a
        PSF with the pixel in which case 'no_pixel' is the right choice.
    src_inds : array-like
        A list of identifiers for each object. The entries from this list are
        fed to `src_func` in order to determine which object to draw.
    src_func : function-like
        A function with signature `gal_func(src_ind)` that returns a tuple of
            src : `galsim.GSObject` to be rendered.
            image_pos : `galsim.PositionD` the position of the object in the
                image
    n_jobs : int, optional
        The number of process to use. If None, then default to the number
        of CPUs as determined by the `loky` package in `joblib.externals`.

    Returns
    -------
    image : galsim.ImageD
        The full image with all of the sources rendered.
    """

    if n_jobs is None:
        n_jobs = joblib.externals.loky.cpu_count()
    n_srcs_per_job = int(np.ceil(len(src_inds) / n_jobs))

    jobs = []
    for job_ind in range(n_jobs):
        start = job_ind * n_srcs_per_job
        end = min(start + n_srcs_per_job, len(src_inds))
        jobs.append(joblib.delayed(_render_list)(
            src_inds[start:end], wcs, draw_method, image_shape, src_func))

    with joblib.Parallel(n_jobs=n_jobs, backend='loky', verbose=0) as p:
        outputs = p(jobs)

    # sum the images
    im = outputs[0]
    if len(outputs) > 1:
        for o in outputs[1:]:
            im += o

    return im


def _render_list(inds, wcs, draw_method, image_shape, src_func):
    im = galsim.ImageD(nrow=image_shape[0], ncol=image_shape[1])
    for ind in inds:
        
        try:
            src, pos = src_func(ind)

            stamp = render_source_in_image(
                source=src,
                local_wcs=wcs.local(image_pos=pos),
                image_pos=pos,
                draw_method=draw_method)

        except galsim.errors.GalSimFFTSizeError:
            print("FFT WAS TOO BIG FOR IND", ind)
            continue

        # intersect and add to total image
        overlap = stamp.bounds & im.bounds
        if overlap.area() > 0:
            im[overlap] += stamp[overlap]

    return im


def render_source_in_image(*, source, image_pos, local_wcs, draw_method):
    """Render a source in a stamp in a larger image.

    Parameters
    ----------
    source : galsim.GSObject
        The source to render. It must have the `drawImage` method.
    image_pos : galsim.PositionD
        The center of the source in the image.
    local_wcs : galsim.LocalWCS
        A local WCS instance to use.
    draw_method : str
        The method used to draw the image. See the docs of `GSObject.drawImage`
        for details and options. Usually 'auto' is correct unless using a
        PSF with the pixel in which case 'no_pixel' is the right choice.

    Returns
    -------
    stamp : galsim.ImageD
        The rendered object in the stamp.
    """
    # pre-draw to get size
    _im = source.drawImage(
        wcs=local_wcs,
        method=draw_method,
        setup_only=True).array
    assert _im.shape[0] == _im.shape[1]

    # lower-left corner
    # the extact math here doesn't matter so much
    # the offset computation takes care of this relative to any x_ll, y_ll
    # we only need to make sure the full object fits on the image
    x_ll = int(image_pos.x - (_im.shape[1] - 1)/2)
    y_ll = int(image_pos.y - (_im.shape[0] - 1)/2)

    # get the offset of the center
    # this is the offset of the image center from the object center
    # galsim renders objects at the image center, so we have to add this
    # offset when rendering
    dx = image_pos.x - (x_ll + (_im.shape[1] - 1)/2)
    dy = image_pos.y - (y_ll + (_im.shape[0] - 1)/2)

    # draw for real
    stamp = source.drawImage(
        nx=_im.shape[1],
        ny=_im.shape[0],
        wcs=local_wcs,
        method=draw_method,
        offset=galsim.PositionD(x=dx, y=dy))
    stamp.setOrigin(x_ll, y_ll)

    return stamp
