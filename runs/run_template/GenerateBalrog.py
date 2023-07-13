from astropy.io import fits
import numpy as np
import os, sys
from tqdm import tqdm
import glob, h5py
import yaml
from sklearn.neighbors import BallTree

name     = os.path.basename(os.path.dirname(__file__))
MCAL_DIR = os.environ['MCAL_DIR']
PATH     = MCAL_DIR + '/' + name
print('GETTING MCAL FILES FROM:')
print(PATH)


tilenames = glob.glob(PATH + '/metacal_*_gplus.fits')
tilenames = [t[-23:-11] for t in tilenames]

print('--------------------------')
print('USING TILES:')
print('--------------------------')
for t in tilenames:
    print(t)
print('--------------------------')

N = len(tilenames)

#Remove tiles without results
tmp = []
for t in tilenames:
    try:
        fits.open(PATH + '/metacal_%s_gplus.fits'%(t))
        tmp.append(t)
    except:
        print("SKIPPING", t)
        
tilenames = tmp


with open(os.path.dirname(__file__) + '/config_plus.yaml', 'r') as fp:
    config = yaml.load(fp, Loader=yaml.Loader)
        

def _get_cat_path(tile, band, mode = 'plus'):
    
    args = {'name' : name,
            'mode' : mode,
            'tile' : tile,
            'band' : band}
    
    
    return os.environ['MCAL_DIR'] + r"/%(name)s/SrcExtractor_%(tile)s_g%(mode)s_%(band)s-cat.fits" % args
    
    
FINAL_CAT = []


Input_catalog = fits.open(os.environ['CATCOSMOS_DIR'] + '/input_cosmos_v4.fits')[1].data
Mask = ((Input_catalog['bdf_hlr'] > config['gal_kws']['size_min']) & 
        (Input_catalog['bdf_hlr'] < config['gal_kws']['size_max']) &
        (Input_catalog['mag_i'] > config['gal_kws']['mag_min']) & 
        (Input_catalog['mag_i'] < config['gal_kws']['mag_max']))

print(np.average(Mask), np.sum(Mask))

Input_catalog = Input_catalog[Mask]

#Loop over tiles and match input and output quantities for objects in each tile
for t in tqdm(tilenames, desc = 'Matching Tiles'):
    
    mcal  = fits.open(PATH + '/metacal_%s_gplus.fits'% t)[1].data
    ID    = mcal['id']
    
    Truth = fits.open(PATH + '/Input_%s_gplus-cat.fits'% t )[1].data
    
    assert 'truedet' not in config['gal_kws']['truth_type'], "Cannot compute survey tranfer func if true detection is used"
    
    SE_plus_r  = fits.open(_get_cat_path(t, 'r', mode = 'plus'))[1].data
    RA, DEC    = SE_plus_r['ALPHA_J2000'], SE_plus_r['DELTA_J2000']
    number     = SE_plus_r['NUMBER']
    
    inds = np.intersect1d(number, ID, return_indices = True)[1]
    
    tree = BallTree(np.vstack([DEC[inds], RA[inds]]).T * np.pi/180, leaf_size=2, metric="haversine")
    d, j = tree.query(np.vstack([Truth['dec'],  Truth['ra']]).T * np.pi/180)

    d, j = d[:, 0], j[:, 0]
    d = d * 180/np.pi * 60*60 #convert to arcsec
    
    Mask = d < 0.5
    j = j[Mask] #Keep only ids below 0.5 arcsec

    Nobj = len(Mask)
        
        
    dtype  = np.dtype([('Truth_ind','>f4'), ('ra', '>f4'),('dec', '>f4'),('true_ra', '>f4'), ('true_dec', '>f4'), 
                       ('true_FLUX_r','>f4'),('true_FLUX_i','>f4'),('true_FLUX_z','>f4'), 
                       ('FLUX_r','>f4'),     ('FLUX_i','>f4'),     ('FLUX_z','>f4'), 
                       ('FLUX_r_ERR','>f4'), ('FLUX_i_ERR','>f4'), ('FLUX_z_ERR','>f4'),
                       ('photoz','>f4'), ('d_arcsec','>f4'), ('detected', '?')])
    
    output = np.zeros(Nobj, dtype = dtype)
    
    output['Truth_ind'] = Truth['ind']
    output['true_ra']   = Truth['ra']
    output['true_dec']  = Truth['dec']
    output['d_arcsec']  = d
    
    output['photoz']    = Input_catalog['photoz'][Truth['ind']]
    output['detected']  = Mask
    
    output['true_FLUX_r'] = Input_catalog['flux_r'][Truth['ind']]
    output['true_FLUX_i'] = Input_catalog['flux_i'][Truth['ind']]
    output['true_FLUX_z'] = Input_catalog['flux_z'][Truth['ind']]
    
    output['ra'][Mask]  = RA[inds][j]
    output['dec'][Mask] = DEC[inds][j]
    
    
    
    output['FLUX_r'][Mask] = mcal['mcal_flux_noshear'][j, 0]
    output['FLUX_i'][Mask] = mcal['mcal_flux_noshear'][j, 1]
    output['FLUX_z'][Mask] = mcal['mcal_flux_noshear'][j, 2]
    
    output['FLUX_r_ERR'][Mask] = np.sqrt(mcal['mcal_flux_cov_noshear'][j, 0, 0])
    output['FLUX_i_ERR'][Mask] = np.sqrt(mcal['mcal_flux_cov_noshear'][j, 1, 1])
    output['FLUX_z_ERR'][Mask] = np.sqrt(mcal['mcal_flux_cov_noshear'][j, 2, 2])
    
    
    output['ra'][np.invert(Mask)]  = np.NaN
    output['dec'][np.invert(Mask)] = np.NaN
    
    output['FLUX_r'][np.invert(Mask)] = np.NaN
    output['FLUX_i'][np.invert(Mask)] = np.NaN
    output['FLUX_z'][np.invert(Mask)] = np.NaN
    
    output['FLUX_r_ERR'][np.invert(Mask)] = np.NaN
    output['FLUX_i_ERR'][np.invert(Mask)] = np.NaN
    output['FLUX_z_ERR'][np.invert(Mask)] = np.NaN
        
    FINAL_CAT.append(output)

    
FINAL_CAT = np.concatenate(FINAL_CAT, axis = 0)


with h5py.File(PATH + '/BalrogOfTheDECADE_Catalog.hdf5', 'w') as f:
    
    for i in FINAL_CAT.dtype.names:
        
        f.create_dataset(i, data = FINAL_CAT[i])