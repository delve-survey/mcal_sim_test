import numpy as np
import galsim


class GaussPixPSF(object):
    """A pixelized Gaussian PSF.

    Parameters
    ----------
    gstd : float, optional
        The standard deviation of the Gaussian distribution for drawing the
        PSF shape.
    fwhm_frac_std : float, optional
        The fractional range allowed for deviations in the PSF FWHM.
    s2n : float, optional
        If not `None`, this option forces the code to add noise to the PSF
        image so that it has total S/N `s2n`.

    Methods
    -------
    getPSF(image_pos)
        Get the PSF represented as an interpolated image at a point.
    """
    def __init__(self, *, gstd=0.01, fwhm_frac_std=0.1, s2n=None):
        self.gstd = gstd
        self.fwhm_frac_std = fwhm_frac_std
        self.s2n = s2n

    def getPSF(self, image_pos, wcs):
        """Get the PSF as an InterpolatedImage

        Parameters
        ----------
        image_pos : galsim.PositionD
            The image position at which to draw the PSF model.
        wcs : galsim.BaseWCS or subclass
            The WCS to use to draw the PSF.

        Returns
        -------
        psf : galsim.InterpolatedImage
            The PSF model.
        """
        wcs = wcs.local(image_pos)

        # we seed with the nearest pixel to make things reproducible
        seed = int(image_pos.x + 0.5) * 4096 + int(image_pos.y + 0.5)
        seed = seed % 2**30
        rng = np.random.RandomState(seed=seed)

        g1 = rng.normal() * self.gstd
        g2 = rng.normal() * self.gstd
        fwhm = (
            rng.uniform(low=-self.fwhm_frac_std, high=self.fwhm_frac_std) +
            1.0) * 0.9
        psf = galsim.Gaussian(fwhm=fwhm).shear(g1=g1, g2=g2).withFlux(1.0)
        psf_im = psf.drawImage(
            nx=69, ny=69, scale=0.125, method='no_pixel').array

        if self.s2n is not None:
            noise_std = np.sqrt(np.sum(psf_im**2)/self.s2n**2)
            psf_im += (rng.normal(size=psf_im.shape) * noise_std)

        psf = galsim.InterpolatedImage(
            galsim.ImageD(psf_im),
            scale=0.125,
            ).withFlux(1.0)
        return psf
