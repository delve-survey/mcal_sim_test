import numpy as np
import ngmix
import galsim
import joblib

import sys
sys.path.insert(0,'/home/secco/delve_cs_test/code/newish_metacal/tests')
sys.path.insert(0,'/Users/secco/Documents/projects/delve_shear/delve_cs_test/code/newish_metacal')
sys.path.insert(0,'/Users/secco/Documents/projects/delve_shear/delve_cs_test/code/newish_metacal/metacal')

#sys.path.insert(0,'/home/dhayaa/Desktop/DECADE/delve_cs_test/code/newish_metacal/tests')
#sys.path.insert(0,'/home/dhayaa/Desktop/DECADE/delve_cs_test/code/newish_metacal')
#sys.path.insert(0,'/home/dhayaa/Desktop/DECADE/delve_cs_test/code/newish_metacal/metacal')

from metacal import MetacalFitter


CONFIG = {
    'metacal': {
        # check for an edge hit
        'bmask_flags': 2**30,

        'metacal_pars': {
            'psf': 'fitgauss',
            'types': ['noshear', '1p', '1m', '2p', '2m'],
        },

        'model': 'gauss',

        'max_pars': {
            'ntry': 2,
            'pars': {
                'method': 'lm',
                'lm_pars': {
                    'maxfev': 2000,
                    'xtol': 5.0e-5,
                    'ftol': 5.0e-5,
                }
            }
        },

        'priors': {
            'cen': {
                'type': 'normal2d',
                'sigma': 0.263
            },

            'g': {
                'type': 'ba',
                'sigma': 0.2
            },

            'T': {
                'type': 'two-sided-erf',
                'pars': [-1.0, 0.1, 1.0e+06, 1.0e+05]
            },

            'flux': {
                'type': 'two-sided-erf',
                'pars': [-100.0, 1.0, 1.0e+09, 1.0e+08]
            }
        },

        'psf': {
            'model': 'gauss',
            'ntry': 2,
            'lm_pars': {
                'maxfev': 2000,
                'ftol': 1.0e-5,
                'xtol': 1.0e-5
            }
        }
    },
}


SHEARS = ['noshear', '1p', '1m', '2p', '2m']


def make_sim(*, seed, g1, g2, s2n=1e6):
    rng = np.random.RandomState(seed=seed)

    gal = galsim.Exponential(half_light_radius=0.5).shear(g1=g1, g2=g2)
    psf = galsim.Gaussian(fwhm=0.9)
    obj = galsim.Convolve([gal, psf])
    dim = 53
    cen = (dim-1)/2
    dither = rng.uniform(size=2, low=-0.5, high=0.5)
    scale = 0.263

    im = obj.drawImage(nx=53, ny=53, offset=dither, scale=scale).array
    nse = np.sqrt(np.sum(im**2)) / s2n
    im += rng.normal(size=im.shape, scale=nse)

    psf_im = psf.drawImage(nx=53, ny=53, scale=scale).array

    jac = ngmix.DiagonalJacobian(
        scale=scale, row=cen+dither[0]+1, col=cen+dither[1]+1
    )
    psf_jac = ngmix.DiagonalJacobian(
        scale=scale, row=cen+1, col=cen+1
    )

    obs = ngmix.Observation(
        image=im,
        weight=np.ones_like(im) / nse**2,
        jacobian=jac,
        bmask=np.zeros_like(im, dtype=np.int32),
        psf=ngmix.Observation(
            image=psf_im,
            jacobian=psf_jac,
        ),
        meta={"orig_row": cen, "orig_col": cen},
    )
    return ngmix.observation.get_mb_obs(obs)


def run_single_sim_pair(seed, s2n=1e6):
    mbobs_plus = make_sim(seed=seed, g1=0.02, g2=0.0, s2n=s2n)
    rng = np.random.RandomState(seed=seed)
    ftr = MetacalFitter(CONFIG, 1, rng)
    ftr.go([mbobs_plus])
    res_p = ftr.result
    if res_p is None:
        return None

    mbobs_minus = make_sim(seed=seed, g1=-0.02, g2=0.0, s2n=s2n)
    rng = np.random.RandomState(seed=seed)
    ftr = MetacalFitter(CONFIG, 1, rng)
    ftr.go([mbobs_minus])
    res_m = ftr.result
    if res_m is None:
        return None

    return res_p, res_m


def _msk_it(*, d, s2n_cut, size_cut, shear=''):
    return (
        (d['mcal_flags'] == 0) &
        (d['mcal_s2n' + shear] > s2n_cut) &
        (d['mcal_T_ratio' + shear] > size_cut)
    )


def measure_g1g2R(*, d, s2n_cut, size_cut):
    msks = {}
    for shear in SHEARS:
        msks[shear] = _msk_it(
            d=d, s2n_cut=s2n_cut, size_cut=size_cut, shear='_' + shear)

    g1_1p = np.mean(d['mcal_g_1p'][msks['1p'], 0])
    g1_1m = np.mean(d['mcal_g_1m'][msks['1m'], 0])
    g2_2p = np.mean(d['mcal_g_2p'][msks['2p'], 1])
    g2_2m = np.mean(d['mcal_g_2m'][msks['2m'], 1])
    R11 = (g1_1p - g1_1m) / 2 / 0.01
    R22 = (g2_2p - g2_2m) / 2 / 0.01

    g1 = np.mean(d['mcal_g_noshear'][msks['noshear'], 0])
    g2 = np.mean(d['mcal_g_noshear'][msks['noshear'], 1])

    return g1, g2, R11, R22


def measure_m_c(res_p, res_m):
    g1p, g2p, R11p, R22p = measure_g1g2R(d=res_p, s2n_cut=10, size_cut=0.5)
    g1m, g2m, R11m, R22m = measure_g1g2R(d=res_m, s2n_cut=10, size_cut=0.5)

    m = (g1p - g1m)/(R11p + R11m)/0.02 - 1
    c = (g2p + g2m)/(R22p + R22m)
    return m, c


def measure_m_c_bootstrap(res_p, res_m, seed, nboot=100):
    rng = np.random.RandomState(seed=seed)
    marr = []
    carr = []
    for _ in range(nboot):
        inds = rng.choice(res_p.shape[0], size=res_p.shape[0], replace=True)
        m, c = measure_m_c(res_p[inds], res_m[inds])
        marr.append(m)
        carr.append(c)

    m, c = measure_m_c(res_p, res_m)
    return m, np.std(marr), c, np.std(carr)


def test_metacal():
    nsims = 100

    rng = np.random.RandomState(seed=34132)
    seeds = rng.randint(size=nsims, low=1, high=2**29)
    jobs = [
        joblib.delayed(run_single_sim_pair)(seed)
        for seed in seeds
    ]
    outputs = joblib.Parallel(n_jobs=-1, verbose=10)(jobs)
    res_p = []
    res_m = []
    for res in outputs:
        if res is not None:
            res_p.append(res[0])
            res_m.append(res[1])
    res_p = np.concatenate(res_p)
    res_m = np.concatenate(res_m)

    seed = rng.randint(size=nsims, low=1, high=2**29)
    m, merr, c, cerr = measure_m_c_bootstrap(res_p, res_m, seed, nboot=100)

    print("m: %f +/- %f [1e-3, 3-sigma]" % (m/1e-3, 3*merr/1e-3), flush=True)
    print("c: %f +/- %f [1e-5, 3-sigma]" % (c/1e-5, 3*cerr/1e-5), flush=True)

    assert np.abs(m) < max(5e-4, 3*merr), (m, merr)
    assert np.abs(c) < 4.0*cerr, (c, cerr)


def test_metacal_slow():
    nsims = 10_000

    rng = np.random.RandomState(seed=342)
    seeds = rng.randint(size=nsims, low=1, high=2**29)
    jobs = [
        joblib.delayed(run_single_sim_pair)(seed, s2n=20)
        for seed in seeds
    ]
    outputs = joblib.Parallel(n_jobs=-1, verbose=10)(jobs)
    res_p = []
    res_m = []
    for res in outputs:
        if res is not None:
            res_p.append(res[0])
            res_m.append(res[1])
    res_p = np.concatenate(res_p)
    res_m = np.concatenate(res_m)

    seed = rng.randint(size=nsims, low=1, high=2**29)
    m, merr, c, cerr = measure_m_c_bootstrap(res_p, res_m, seed, nboot=100)

    print("m: %f +/- %f [1e-3, 3-sigma]" % (m/1e-3, 3*merr/1e-3), flush=True)
    print("c: %f +/- %f [1e-5, 3-sigma]" % (c/1e-5, 3*cerr/1e-5), flush=True)

    assert np.abs(m) < max(5e-4, 3*merr), (m, merr)
    assert np.abs(c) < 4.0*cerr, (c, cerr)


if __name__ == "__main__":

	test_metacal_slow()
