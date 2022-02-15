import numpy as np
import galsim


class NonGaussPixPSF(object):
    """A pixelized non-Gaussian PSF (galsim's Moffat profile).
    Moffat's surface brightness profile goes as I(r)~(1+(r/r_0)^2)^beta

    Parameters
    ----------
    gstd : float, optional
        A small shear applied to the PSF shape to be drawn
    beta : float
        The beta parameter in the Moffat profile
    scale_radius : float
        The scale radius  distribution for drawing the Moffat PSF shape (in arcsec, hopefully).
    fwhm_frac_std : float, optional
        The fractional range allowed for deviations in the PSF scale radius.
    s2n : float, optional
        If not `None`, this option forces the code to add noise to the PSF
        image so that it has total S/N `s2n`.

    Methods
    -------
    getPSF(image_pos)
        Get the PSF represented as an interpolated image at a point.
    """
    def __init__(self, *,draw_with_wcs,gstd=0.01, beta=3.0, scale_radius=0.9,fwhm_frac_std=0.1, s2n=None):
        self.gstd = gstd
        self.beta = beta
        self.scale_radius = scale_radius
        self.fwhm_frac_std = fwhm_frac_std
        self.s2n = s2n
        self.draw_with_wcs = draw_with_wcs
        print('\nFrom nongauss_pix_psf/NonGaussPixPSF: draw_with_wcs=',draw_with_wcs,'\n')

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
            1.0) * self.scale_radius
        #psf = galsim.Gaussian(fwhm=fwhm).shear(g1=g1, g2=g2).withFlux(1.0)
        psf = galsim.Moffat(beta=self.beta,scale_radius=fwhm).shear(g1=g1, g2=g2).withFlux(1.0)
        
        if self.draw_with_wcs==False:
            psf_im = psf.drawImage(
                nx=69, ny=69, scale=0.125, method='no_pixel').array
        if self.draw_with_wcs==True:
             psf_im = psf.drawImage(
                nx=69, ny=69, wcs=wcs, method='no_pixel').array

        if self.s2n is not None:
            noise_std = np.sqrt(np.sum(psf_im**2)/self.s2n**2)
            psf_im += (rng.normal(size=psf_im.shape) * noise_std)

        if self.draw_with_wcs==False:
            psf = galsim.InterpolatedImage(
                galsim.ImageD(psf_im),
                scale=0.125,
                ).withFlux(1.0)
        if self.draw_with_wcs==True:
            psf = galsim.InterpolatedImage(
                galsim.ImageD(psf_im),
                wcs=wcs,
                ).withFlux(1.0)   

        return psf
