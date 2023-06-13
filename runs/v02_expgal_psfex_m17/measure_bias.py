from astropy.io import fits
import numpy as np
import os
from tqdm import tqdm
import sys
import yaml

from SimButler import tilenames


summary = np.mean

tilenames = tilenames

print('--------------------------')
print('USING TILES:')
print('--------------------------')
for t in tilenames:
    print(t)
print('--------------------------')

N = len(tilenames)

name     = os.path.basename(os.path.dirname(__file__))
MCAL_DIR = os.environ['MCAL_DIR']
PATH     = MCAL_DIR + '/' + name
print('GETTING MCAL FILES FROM:')
print(PATH)

tiles_list = [[x] for x in tilenames]
tiles_list.append(tilenames)

seed_list = [[i] for i in range(len(tilenames))]
seed_list.append([i for i in range(len(tilenames))])

tiles_list = [tiles_list[-1]]
seed_list  = [seed_list[-1]]

with open(os.path.dirname(__file__) + '/config_plus.yaml', 'r') as fp:
    config = yaml.load(fp, Loader=yaml.Loader)
        
def _get_subset_array(res, i, string):
    '''
    Convenience function for getting reading out specific arrays
    and only for specific objects that are part of a jackknife "patch"
    '''
    ind = np.where(Results['batch_num_%s'%string[3:]] != i)[0] #Select only objects in given batch/subset
    
    return res[string][ind]

def _get_cat_path(tile, band, seed, mode = 'plus'):
    
    args = {'name' : name,
            'seed' : seed,
            'mode' : mode,
            'tile' : tile,
            'band' : band}
    
    
    return os.environ['MCAL_DIR'] + r"/%(name)s/SrcExtractor_%(tile)s_g%(mode)s_%(band)s-cat.fits" % args
    
#Loop over tiles (and also a combination of all tiles together)
#and get m and c for objects in each tile
for tiles, seed in zip(tiles_list, seed_list):
    
    _t = []
    _s = []
    
    for t, i in zip(tiles, seed):
        
        try:
            fits.open(PATH + '/metacal_%s_gplus.fits'%(t))
            fits.open(PATH + '/metacal_%s_gminus.fits'%(t))
            _t.append(t)
            _s.append(i)
        except:
            print("SKIPPING", t)
        
        

    tiles = _t
    seed  = _s

    if len(tiles) == 0: continue

    N = len(tiles) 
    print("N_TILES:", N)    
    #if N == 1: continue #Just quick "hack" so we skip all single tile calculations
    
    print("--------------------------------------------------------")
    print("TILES:", tiles)
    print("--------------------------------------------------------")
    
    gplus  = [fits.open(PATH + '/metacal_%s_gplus.fits'%(t))[1].data for t, i in zip(tiles, seed)]
    gminus = [fits.open(PATH + '/metacal_%s_gminus.fits'%(t))[1].data for t, i in zip(tiles, seed)]
    
    if 'truedet' not in config['gal_kws']['truth_type']:
        print(_get_cat_path(tiles[0], 'r', seed[0], mode = 'plus'))
        SE_plus_r  = [fits.open(_get_cat_path(t, 'r', i, mode = 'plus'))[1].data for t, i in zip(tiles, seed)]
        SE_plus_i  = [fits.open(_get_cat_path(t, 'i', i, mode = 'plus'))[1].data for t, i in zip(tiles, seed)]
        SE_plus_z  = [fits.open(_get_cat_path(t, 'z', i, mode = 'plus'))[1].data for t, i in zip(tiles, seed)]

        SE_minus_r = [fits.open(_get_cat_path(t, 'r', i, mode = 'minus'))[1].data for t, i in zip(tiles, seed)]
        SE_minus_i = [fits.open(_get_cat_path(t, 'i', i, mode = 'minus'))[1].data for t, i in zip(tiles, seed)]
        SE_minus_z = [fits.open(_get_cat_path(t, 'z', i, mode = 'minus'))[1].data for t, i in zip(tiles, seed)]
        
        
        mask_SEflag_plus = ((np.concatenate([SE_plus_r[i]['FLAGS'] for i in range(N)]) <= 3) &
                            (np.concatenate([SE_plus_i[i]['FLAGS'] for i in range(N)]) <= 3) &
                            (np.concatenate([SE_plus_z[i]['FLAGS'] for i in range(N)]) <= 3))
        
        mask_snr_plus = (np.concatenate([SE_plus_r[i]['FLUX_AUTO'] for i in range(N)])/
                         np.concatenate([SE_plus_r[i]['FLUXERR_AUTO'] for i in range(N)])) > 5
        
        mask_size_plus = (np.concatenate([SE_plus_r[i]['FLUX_RADIUS']*0.236 for i in range(N)]) > 0.5)
        
        
        mask_SEflag_minus = ((np.concatenate([SE_minus_r[i]['FLAGS'] for i in range(N)]) <= 3) &
                             (np.concatenate([SE_minus_i[i]['FLAGS'] for i in range(N)]) <= 3) &
                             (np.concatenate([SE_minus_z[i]['FLAGS'] for i in range(N)]) <= 3))
        
        mask_snr_minus = (np.concatenate([SE_minus_r[i]['FLUX_AUTO'] for i in range(N)])/
                          np.concatenate([SE_minus_r[i]['FLUXERR_AUTO'] for i in range(N)])) > 5
        
        mask_size_minus = (np.concatenate([SE_minus_r[i]['FLUX_RADIUS']*0.236 for i in range(N)]) > 0.5)

        
        number_plus  = np.concatenate([SE_plus_r[i]['NUMBER'] + i*int(1e6) for i in range(N)])
        number_minus = np.concatenate([SE_minus_r[i]['NUMBER'] + i*int(1e6) for i in range(N)])
        print("Size before cuts", len(number_plus))
        
    ngrid = 10
    
    spacing = 10_000/ngrid
    Results = {}
    
    for g_name, g in zip(['p', 'm'], [gplus, gminus]):
        
        Flags = np.concatenate([g[i]['mcal_flags'] for i in range(N)])
        ID    = np.concatenate([g[i]['id'] + i*int(1e6) for i in range(N)])
        bfrac = np.concatenate([g[i]['badfrac'] for i in range(N)])

        for s in ['noshear', '1p', '1m', '2p', '2m']:
              
            
            SNR    = np.concatenate([g[i]['mcal_s2n_%s'%s] for i in range(N)])
            Tratio = np.concatenate([g[i]['mcal_T_ratio_%s'%s] for i in range(N)])
            g1, g2 = np.concatenate([g[i]['mcal_g_%s'%s] for i in range(N)]).T
            T      = np.concatenate([g[i]['mcal_T_%s'%s] for i in range(N)])
            mag_r  = 30 - 2.5*np.log10(np.concatenate([g[i]['mcal_flux_%s'%s][:, 0] for i in range(N)]))

            Mask = (Flags == 0) & (SNR > 10) & (SNR < 1e10) & (T < 10) & (Tratio > 0.5) #& (bfrac < 1) & np.invert((T > 2) & (SNR < 30)) & np.invert((np.log10(T) < (22.25 - 0)/3.5) & (g1**2 + g2**2 > 0.8**2))
             
            if 'truedet' not in config['gal_kws']['truth_type']:
                
                if g_name == 'p':
                    inds = np.intersect1d(ID, number_plus, return_indices = True)[2]
                    Mask = Mask & (mask_SEflag_plus & mask_snr_plus & mask_size_plus)[inds]
                elif g_name == 'm':
                    inds = np.intersect1d(ID, number_minus, return_indices = True)[2]
                    Mask = Mask & (mask_SEflag_minus & mask_snr_minus & mask_size_minus)[inds]
            
            Results['e1_%s_%s'%(g_name, s)] = np.concatenate([g[i]['mcal_g_%s'%s][:,0] for i in range(N)])[Mask]
            Results['e2_%s_%s'%(g_name, s)] = np.concatenate([g[i]['mcal_g_%s'%s][:,1] for i in range(N)])[Mask]
            
            Results['batch_num_%s_%s'%(g_name, s)] = np.concatenate([
                                                        g[i]['x']//spacing + 
                                                        g[i]['y']//spacing * ngrid + 
                                                        i * ngrid**2 
                                                        for i in range(N)]).astype(int)[Mask]

    del gplus, gminus
    

    Npatch = ngrid**2*len(tilenames)
    N = np.arange(Npatch)
#     N = np.random.choice(N, 200, replace = False)
    e1_plus  = np.zeros(N.size)
    e1_minus = np.zeros(N.size)
    e2_plus  = np.zeros(N.size)
    e2_minus = np.zeros(N.size)

    R11_plus  = np.zeros(N.size)
    R11_minus = np.zeros(N.size)
    R22_plus  = np.zeros(N.size)
    R22_minus = np.zeros(N.size)

    
    for i in tqdm(range(N.size)):

        e1_plus[i]  = summary(_get_subset_array(Results, N[i], 'e1_p_noshear'))
        e1_minus[i] = summary(_get_subset_array(Results, N[i], 'e1_m_noshear'))
        e2_plus[i]  = summary(_get_subset_array(Results, N[i], 'e2_p_noshear'))
        e2_minus[i] = summary(_get_subset_array(Results, N[i], 'e2_m_noshear'))
        
        R11_plus[i]  = (summary(_get_subset_array(Results, N[i], 'e1_p_1p')) - summary(_get_subset_array(Results, N[i], 'e1_p_1m')))/0.02
        R11_minus[i] = (summary(_get_subset_array(Results, N[i], 'e1_m_1p')) - summary(_get_subset_array(Results, N[i], 'e1_m_1m')))/0.02
        R22_plus[i]  = (summary(_get_subset_array(Results, N[i], 'e2_p_2p')) - summary(_get_subset_array(Results, N[i], 'e2_p_2m')))/0.02
        R22_minus[i] = (summary(_get_subset_array(Results, N[i], 'e2_m_2p')) - summary(_get_subset_array(Results, N[i], 'e2_m_2m')))/0.02

    R11 = (R11_plus + R11_minus)/2
    R22 = (R22_plus + R22_minus)/2

    g1 = (e1_plus - e1_minus)/(2*R11)
    g2 = (e2_plus + e2_minus)/(2*R22)

    m = g1/0.02 - 1
    c = g2

    print('\n---------------------------------\n')
    print("----------- COMBINED ---------------")
    print('size      N = %d (Masked %0.2f)'%(Results['e1_p_noshear'].shape[0], np.average(Mask)))
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))
    
    print('\n---------------------------------\n')
    
    with open(os.path.dirname(__file__) + '/README.md', 'w') as f:
        
        f.write("TEST RESULTS FOR %s\n" %name)
        f.write("---------------------------------\n")
        f.write("size      N = %d (Masked %0.2f)\n"%(Results['e1_p_noshear'].shape[0], np.average(Mask)))
        f.write("recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]\n"%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
        f.write("recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]\n"%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))
        
    R11 = R11_plus
    R22 = R22_plus

    g1 = e1_plus/R11
    g2 = e2_plus/R22
    
    m = g1/0.02 - 1
    c = g2

    print("----------- PLUS ---------------")
    print('size      N = %d'%Results['e1_p_noshear'].shape[0])
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))
    
    R11 = R11_minus
    R22 = R22_minus

    g1 = e1_minus/R11
    g2 = e2_minus/R22
    
    m = g1/(-0.02) - 1
    c = g2

    print("----------- MINUS ---------------")
    print('size      N = %d'%Results['e1_p_noshear'].shape[0])
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))
          
          
    print("------------------------------------------------------------")
