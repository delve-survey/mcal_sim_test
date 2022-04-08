import os

import numpy as np
import galsim
from galsim.des import DES_PSFEx

class PSFEx_Deconv(object):
    """A wrapper for PSFEx to use with Galsim. Main use is deconvolving pixel scale.

    This wrapper just takes in galsim PSFEx and deconvolves it with pixel scale

    Parameters
    ----------
    file_name : str
        The file with the Psfex psf solution.
    """
    
    _req_params = {'file_name': str}
    _opt_params = {}
    _single_params = []
    _takes_rng = False

    def __init__(self, file_name, wcs):
        self.file_name = file_name
        self.wcs = wcs
        self._psfex = DES_PSFEx(os.path.expanduser(os.path.expandvars(file_name)), wcs = wcs)

    def getPSFEx(self):
        return self._psfex
    
    def getPSF(self, image_pos, wcs=None):
        """Get a deconvolved image of the PSF at the given location.

        Parameters
        ----------
        image_pos : galsim.Position
            The image position for the PSF.
        wcs : galsim.BaseWCS or subclass
            The WCS to use to draw the PSF.

        Returns
        -------
        psf : galsim.InterpolatedImage
            The PSF at the image position.
        """
        
        x_interpolant='lanczos15'
        gsparams=None
        
        scale = 0.25
        pixel_wcs = galsim.PixelScale(scale)

        # nice and big image size here cause this has been a problem
#         image = galsim.ImageD(ncol=19, nrow=19, wcs=pixel_wcs)

#         psf = self.getPSFEx().getPSF(image_pos).drawImage(
#             #center=image_pos, #Dhayaa: Not using center here to drae image. PSF is already obtained at image_pos.
#             image=image,
#             offset=(0, 0)) #explicitly setting to zero because we don't need offset

#         psf = galsim.InterpolatedImage(
#             galsim.ImageD(psf.array),  # make sure galsim is not keeping state
#             wcs=pixel_wcs,
#             gsparams=gsparams,
#             x_interpolant=x_interpolant
#         )

        psf = self.getPSFEx().getPSF(image_pos)

        psf = galsim.Convolve(
            [psf, galsim.Deconvolve(galsim.Pixel(scale))]
        ).withFlux(
            1.0
        )

        return psf