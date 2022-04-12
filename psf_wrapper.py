import galsim
import galsim.des

#from ..des_piff import DES_Piff
from gauss_pix_psf import GaussPixPSF
from nongauss_pix_psf import NonGaussPixPSF
from galsim.des import DES_PSFEx
from des_psfex import DES_PSFEx_Deconv
from psfex_deconvolved import PSFEx_Deconv

class PSFWrapper(object):
    """Wrapper to interface galsim objects.

    This class combines APIs from the `psfex` Python package (used by the
    MEDS code) and the galsim.des.DES_PSFEx module (used for image
    simulations).

    Parameters
    ----------
    psf : galsim.GSObject or a subclass
        The galsim object to draw.
    wcs : a galsim WCS object
        The WCS to use to convert the PSF from world coordinates to image
        coordinates.
    n_pix : int
        The number of pixels on a side for PSF image. Make sure to make it
        an odd number.

    Methods
    -------
    get_rec(row, col)
        Get a reconstruction of the PSF.
    get_rec_shape(row, col)
        Get the image shape of a reconstruction of the PSF.
    getPSF(image_pos)
        Get the PSF as a galsim.GSObject at a given image position.
    get_center(row, col)
        Get the center of the PSF in the stamp/cutout.
    get_sigma(row, col)
        Raises a `NotImplementedError`. Here to make sure no code is using it.
    """
    def __init__(self, psf, wcs, n_pix=53):
        self.psf = psf
        self.wcs = wcs
        self.n_pix = n_pix

    def get_rec_shape(self, row, col):
        """Get the shape of the PSF image at a position.

        Parameters
        ----------
        row : float
            The row at which to get the PSF image in the stamp in
            zero-offset image coordinates.
        col : float
            The col at which to get the PSF image in the stamp in
            zero-offset image coordinates.

        Returns
        -------
        psf_shape : tuple of ints
            The shape of the PSF image.
        """
        return (self.n_pix, self.n_pix)

    def getPSF(self, image_pos):
        """Get the PSF as a galsim.GSObject at a given image position.

        Parameters
        ----------
        image_pos : galsim.PositionD
            The image position in one-indexed, pixel centered coordinates.

        Returns
        -------
        psf : galsim.GSObject
            The PSF as a galsim.GSOjbect.
        """
        if isinstance(self.psf, galsim.GSObject):
            return self.psf
        
        #elif isinstance(self.psf, DES_Piff):
        #    wcs = self.wcs.local(image_pos)
        #    return self.psf.getPSF(image_pos, wcs)

        elif isinstance(self.psf, NonGaussPixPSF):
            wcs = self.wcs.local(image_pos)
            return self.psf.getPSF(image_pos, wcs)
        
        elif isinstance(self.psf, GaussPixPSF):
            wcs = self.wcs.local(image_pos)
            return self.psf.getPSF(image_pos, wcs)
        
        elif isinstance(self.psf, DES_PSFEx):
            return self.psf.getPSF(image_pos) #Wrapper doesn't take wcs. Need to pass it when reading file.
        
        elif isinstance(self.psf, DES_PSFEx_Deconv):
            wcs = self.wcs.local(image_pos)
            return self.psf.getPSF(image_pos, wcs)
        
        elif isinstance(self.psf, PSFEx_Deconv):
            wcs = self.wcs.local(image_pos)
            return self.psf.getPSF(image_pos, wcs)
        
        else:
            raise ValueError(
                'We did not recognize the PSF type! %s' % self.psf)

    def get_rec(self, row, col):
        """Get the PSF at a position.

        Parameters
        ----------
        row : float
            The row at which to get the PSF image in the stamp in
            zero-offset image coordinates.
        col : float
            The col at which to get the PSF image in the stamp in
            zero-offset image coordinates.

        Returns
        -------
        psf : np.ndarray, shape (npix, npix)
            An image of the PSF.
        """
        # we add 1 to the positions here since the MEDS code uses
        # zero offset positions and galsim + DES stuff expects one-offset
        im_pos = galsim.PositionD(x=col+1, y=row+1)
        wcs = self.wcs.local(im_pos)
        
        if isinstance(self.psf, galsim.GSObject):
            psf_im = self.psf.drawImage(
                nx=self.n_pix, ny=self.n_pix,
                wcs=wcs).array
        
        #elif isinstance(self.psf, DES_Piff):
        #    psf_at_pos = self.psf.getPSF(im_pos, wcs)
        #    psf_im = psf_at_pos.drawImage(
        #        wcs=wcs, nx=self.n_pix, ny=self.n_pix).array
 
        elif isinstance(self.psf, NonGaussPixPSF):
            psf_at_pos = self.psf.getPSF(im_pos, wcs)
            psf_im = psf_at_pos.drawImage(
                wcs=wcs, nx=self.n_pix, ny=self.n_pix).array      

        elif isinstance(self.psf, GaussPixPSF):
            psf_at_pos = self.psf.getPSF(im_pos, wcs)
            psf_im = psf_at_pos.drawImage(
                wcs=wcs, nx=self.n_pix, ny=self.n_pix).array
        
        elif isinstance(self.psf, DES_PSFEx):
            psf_at_pos = self.psf.getPSF(im_pos) #No wcs passed here. Need to pass when reading file.
            psf_im = psf_at_pos.drawImage(
                wcs=wcs, nx=self.n_pix, ny=self.n_pix, method = 'no_pixel').array #Need to use no_pixel because PSF is already convolved with pixel scale
        
        elif isinstance(self.psf, DES_PSFEx_Deconv):
            psf_at_pos = self.psf.getPSF(im_pos) #No wcs passed here. Need to pass when reading file.
            psf_im = psf_at_pos.drawImage(
                wcs=wcs, nx=self.n_pix, ny=self.n_pix).array 
            
        elif isinstance(self.psf, PSFEx_Deconv):
            psf_at_pos = self.psf.getPSF(im_pos, self.wcs) #Need to pass wcs here in order to get local jacobian and transform pixel scale
            psf_im = psf_at_pos.drawImage(
                wcs=wcs, nx=self.n_pix, ny=self.n_pix).array 
        
        else:
            raise ValueError(
                'We did not recognize the PSF type! %s' % self.psf)

        # commented out to make sure this is never done
        # usually this does not help anything
        # leaving notes here for the scientists of the future
        # if self.snr is not None:
        #     npix = psf_im.shape[0] * psf_im.shape[1]
        #     sigma = psf_im.sum() / (self.snr * npix**0.5)
        #     print("adding psf noise, final psf s/n = %f" % (
        #         psf_im.sum()/sigma**2 / np.sqrt(npix/sigma**2)))
        #     noise = np.random.normal(scale=sigma, size=psf_im.shape)
        #     psf_im += noise

        return psf_im

    def get_center(self, row, col):
        """Get the center of the PSF in the stamp/cutout.

        Parameters
        ----------
        row : float
            The row at which to get the PSF center in the stamp in
            zero-offset image coordinates.
        col : float
            The col at which to get the PSF center in the stamp in
            zero-offset image coordinates.

        Returns
        -------
        cen : 2-tuple of floats
            The center of the PSF in zero-offset image coordinates.
        """
        return (self.n_pix-1.)/2., (self.n_pix-1.)/2.

    def get_sigma(self, row, col):
        # note this used to return -99
        raise NotImplementedError()
