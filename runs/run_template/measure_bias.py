from astropy.io import fits
import numpy as np
import os
from tqdm import tqdm
import sys

from render_jobs import tilenames


summary = np.mean

tilenames = tilenames[:5]

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

def _get_subset_array(res, i, string):
    '''
    Convenience function for getting reading out specific arrays
    and only for specific objects that are part of a jackknife "patch"
    '''
    ind = np.where(Results['batch_num_%s'%string[3:]] != i)[0] #Select only objects in given batch/subset
    
    return res[string][ind]

#Loop over tiles (and also a combination of all tiles together)
#and get m and c for objects in each tile
for tiles, seed in zip(tiles_list, seed_list):
    
    N = len(tiles) 
    
    if N == 1: continue #Just quick "hack" so we skip all single tile calculations
    
    print("--------------------------------------------------------")
    print("TILES:", tiles)
    print("--------------------------------------------------------")
    
    gplus  = [fits.open(PATH + '/metacal_%s_seed%d_gplus.fits'%(t, i)) for t, i in zip(tiles, seed)]
    gminus = [fits.open(PATH + '/metacal_%s_seed%d_gminus.fits'%(t, i)) for t, i in zip(tiles, seed)]
    
    
    ngrid = 10
    
    spacing = 10_000/ngrid
    Results = {}
    
    for g_name, g in zip(['p', 'm'], [gplus, gminus]):
        
        Flags = np.concatenate([g[i][1].data['mcal_flags'] for i in range(N)])

        for s in ['noshear', '1p', '1m', '2p', '2m']:
              
            
            SNR = np.concatenate([g[i][1].data['mcal_s2n_%s'%s] for i in range(N)])
            T   = np.concatenate([g[i][1].data['mcal_T_ratio_%s'%s] for i in range(N)])
            Mask = (Flags == 0) & (SNR > 10) & (T > 0.5)
            
            Results['e1_%s_%s'%(g_name, s)] = np.concatenate([g[i][1].data['mcal_g_%s'%s][:,0] for i in range(N)])[Mask]
            Results['e2_%s_%s'%(g_name, s)] = np.concatenate([g[i][1].data['mcal_g_%s'%s][:,1] for i in range(N)])[Mask]
            
            Results['batch_num_%s_%s'%(g_name, s)] = np.concatenate([
                                                        g[i][1].data['x']//spacing + 
                                                        g[i][1].data['y']//spacing * ngrid + 
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
    print('size      N = %d'%Results['e1_p_noshear'].shape)
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))
    
    print('\n---------------------------------\n')
    
    R11 = R11_plus
    R22 = R22_plus

    g1 = e1_plus/R11
    g2 = e2_plus/R22
    
    m = g1/0.02 - 1
    c = g2

    print("----------- PLUS ---------------")
    print('size      N = %d'%Results['e1_p_noshear'].shape)
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))
    
    R11 = R11_minus
    R22 = R22_minus

    g1 = e1_minus/R11
    g2 = e2_minus/R22
    
    m = g1/(-0.02) - 1
    c = g2

    print("----------- MINUS ---------------")
    print('size      N = %d'%Results['e1_p_noshear'].shape)
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))
          
          
    print("------------------------------------------------------------")
