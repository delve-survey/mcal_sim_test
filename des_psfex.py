import os

import numpy as np
import galsim
from galsim.des import DES_PSFEx
import ngmix
from ngmix.fitting import Fitter as LMSimple
from ngmix.admom import AdmomFitter as Admom

from scipy.interpolate import CloughTocher2DInterpolator



class DES_PSFEx_Deconv(object):
    """A wrapper for PSFEx to use with Galsim. Main use is deconvolving pixel scale.

    This wrapper uses ngmix to fit smooth models to the Piff PSF images. The
    parameters of these models are then interpolated across the SE image
    and used to generate a smooth approximation to the PSF.

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
        self._did_fit = False

    def _fit_smooth_model(self):
        dxy = 256
        ny = 4096 // dxy + 1
        nx = 2048 // dxy + 1

        xloc = np.empty((ny, nx), dtype=np.float64)
        yloc = np.empty((ny, nx), dtype=np.float64)
        pars = np.empty((ny, nx, 3), dtype=np.float64)
        for yi, yl in enumerate(np.linspace(1, 4096, ny)):
            for xi, xl in enumerate(np.linspace(1, 2048, nx)):
                rng = np.random.RandomState(seed=yi + nx * xi)
                xloc[yi, xi] = xl
                yloc[yi, xi] = yl

                pos = galsim.PositionD(x=xl, y=yl)
                img = self._draw(pos).drawImage(
                    nx=19, ny=19, scale=0.25, method='sb').array
                nse = np.std(
                    np.concatenate([img[0, :], img[-1, :]]))
                obs = ngmix.Observation(
                    image=img,
                    weight=np.ones_like(img)/nse**2,
                    jacobian=ngmix.jacobian.DiagonalJacobian(
                        x=9, y=9, scale=0.25))

                _g1 = np.nan
                _g2 = np.nan
                _T = np.nan

                
                for count in range(5):
                    try:
                        am = Admom(rng=rng)
                        res = am.go(obs, 0.3)
                        if res['flags'] != 0:
                            continue

                        lm = LMSimple('turb')
                        lm_res = lm.go(obs, res['pars'])
                        if lm_res['flags'] == 0:
                            _g1 = lm_res['pars'][2]
                            _g2 = lm_res['pars'][3]
                            _T = lm_res['pars'][4]
                            break
                    except ngmix.gexceptions.GMixRangeError:
                        pass

                pars[yi, xi, 0] = _g1
                pars[yi, xi, 1] = _g2
                pars[yi, xi, 2] = _T

        xloc = xloc.ravel()
        yloc = yloc.ravel()
        pos = np.stack([xloc, yloc], axis=1)
        assert pos.shape == (xloc.shape[0], 2)

        # make interps
        g1 = pars[:, :, 0].ravel()
        msk = np.isfinite(g1)
        if len(msk) < 10:
            raise ValueError('DES Piff fitting failed too much!')
        if np.any(~msk):
            g1[~msk] = np.mean(g1[msk])
        self._g1int = CloughTocher2DInterpolator(pos, g1)

        g2 = pars[:, :, 1].ravel()
        msk = np.isfinite(g2)
        if len(msk) < 10:
            raise ValueError('DES Piff fitting failed too much!')
        if np.any(~msk):
            g2[~msk] = np.mean(g2[msk])
        self._g2int = CloughTocher2DInterpolator(pos, g2)

        T = pars[:, :, 2].ravel()
        msk = np.isfinite(T)
        if len(msk) < 10:
            raise ValueError('DES Piff fitting failed too much!')
        if np.any(~msk):
            T[~msk] = np.mean(T[msk])
        self._Tint = CloughTocher2DInterpolator(pos, T)

        self._did_fit = True

    def _draw(self, image_pos, x_interpolant='lanczos15', gsparams=None):
        """Get an image of the PSF at the given location.

        Parameters
        ----------
        image_pos : galsim.Position
            The image position for the PSF.
        wcs : galsim.BaseWCS or subclass
            The WCS to use to draw the PSF.
        x_interpolant : str, optional
            The interpolant to use.
        gsparams : galsim.GSParams, optional
            Ootional galsim configuration data to pass along.

        Returns
        -------
        psf : galsim.InterpolatedImage
            The PSF at the image position.
        """
        scale = 0.25
        pixel_wcs = galsim.PixelScale(scale)

        # nice and big image size here cause this has been a problem
        image = galsim.ImageD(ncol=19, nrow=19, wcs=pixel_wcs)

        # piff offsets the center of the PSF from the true image
        # center - here we will return a properly centered image by undoing
        # the offset
        dx = image_pos.x - int(image_pos.x + 0.5)
        dy = image_pos.y - int(image_pos.y + 0.5)

        psf = self.getPSFEx().getPSF(image_pos).drawImage(
            #center=image_pos, #Dhayaa: Not using center here to drae image. PSF is already obtained at image_pos.
            image=image,
            offset=(-dx, -dy))

        psf = galsim.InterpolatedImage(
            galsim.ImageD(psf.array),  # make sure galsim is not keeping state
            wcs=pixel_wcs,
            gsparams=gsparams,
            x_interpolant=x_interpolant
        )

        psf = galsim.Convolve(
            [psf, galsim.Deconvolve(galsim.Pixel(scale))]
        ).withFlux(
            1.0
        )

        return psf

    def getPSFEx(self):
        return self._psfex

    def getPSF(self, image_pos, wcs=None):
        """Get an image of the PSF at the given location.

        Parameters
        ----------
        image_pos : galsim.Position
            The image position for the PSF.
        wcs : galsim.BaseWCS or subclass, optional
            The WCS to use to draw the PSF. Currently ignored.

        Returns
        -------
        psf : galsim.GSObject
            The PSF at the image position.
        """
        if not self._did_fit:
            self._fit_smooth_model()

        arr = np.array([
            np.clip(image_pos.x, 1, 2048),
            np.clip(image_pos.y, 1, 4096)])

        _g1 = self._g1int(arr)[0]
        _g2 = self._g2int(arr)[0]
        _T = self._Tint(arr)[0]
        if np.any(np.isnan(np.array([_g1, _g2, _T]))):
            print("\n\n\n", image_pos, _g1, _g2, _T, "\n\n\n", flush=True)
        pars = np.array([0, 0, _g1, _g2, _T, 1])
        obj = ngmix.gmix.make_gmix_model(pars, 'turb').make_galsim_object()
        return obj.withFlux(1)
