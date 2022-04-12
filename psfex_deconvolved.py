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
        
        psf = self.getPSFEx().getPSF(image_pos) #Get galsim PSF object
        
        wcs = wcs.jacobian(image_pos=image_pos)
        pixel = wcs.toWorld(galsim.Pixel(scale=1)) #Get pixel profile in correct wcs

        deconvolution_kernel = galsim.Deconvolve(pixel) #Create kernel to deconvolve pixel window
        
        psf = galsim.Convolve([psf, deconvolution_kernel]).withFlux(1.0) #Deconvolve

        return psf
